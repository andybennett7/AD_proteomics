options(java.parameters = "-Xmx16g")
# Library
library(tidyverse)
library(ggplot2)
library(readxl)
library(xlsx)
library("org.Hs.eg.db")
library(gplots) #function venn
library(ggvenn)
library(UpSetR)
library(ComplexHeatmap)
library(e1071)
library(dplyr)
library(factoextra)
library(here)
library(pheatmap)
library(pROC)
library(ggpubr)
library(ggsci)
library(pals)
library(gt)
library(ggtext)

getwd()

#### Data Import/Clean/Stats ####
{
  # Proteome Discoverer results
  brain <- read_excel("data/AD_brain_tissue_proteomics_Proteins_W_Isoforms_separator.xlsx", sheet = "Proteins")
  CSF <- read_excel("data/AD_CSF_proteomics_W_Isoforms_separator.xlsx", sheet = "Proteins")
  plasma <- read_excel("data/AD Plasma Proteins Reviewed_WIsoforms_separator.xlsx", sheet = "Proteins")
  
  #depletion list
  high_abundance <- read.csv("data/HI_ABUNDANT_IDS_depletion.csv")
  high_list <- high_abundance$Accession
  
  # Removes columns not generally used
  brain2<- brain %>% dplyr::select(-contains(c("Checked","Exp. q-value", "# PSMs", 
                                               "# AAs", "MW [kDa]", "calc. pI", 
                                               "Score Sequest", "Peptides (by Search Engine):", 
                                               "# Razor Peptides", "# Protein Groups", "Found", 
                                               "Master", "Marked As")))
  CSF2<- CSF %>% dplyr::select(-contains(c("Checked","Exp. q-value", "# PSMs", 
                                           "# AAs", "MW [kDa]", "calc. pI", 
                                           "Score Sequest", "Peptides (by Search Engine):", 
                                           "# Razor Peptides", "# Protein Groups", "Found", 
                                           "Master", "Marked As")))
  plasma2<- plasma %>% dplyr::select(-contains(c("Checked","Exp. q-value", "# PSMs", 
                                                 "# AAs", "MW [kDa]", "calc. pI", 
                                                 "Score Sequest", "Peptides (by Search Engine):", 
                                                 "# Razor Peptides", "# Protein Groups", 
                                                 "Found", "Master", "Marked As")))
  
  
  # smpl start column set to column 8. CHECK THIS!!!
  smpl_start <- 8 
  column_sums_B <- colSums(brain2[smpl_start:ncol(brain2)], na.rm = TRUE)
  
  #smpl_C <- readline(prompt = "Column of first sample?: " )
  column_sums_C <- colSums(CSF2[smpl_start:ncol(CSF2)], na.rm = TRUE)
  
  #smpl_P <- readline(prompt = "Column of first sample?: " )
  column_sums_P <- colSums(plasma2[smpl_start:ncol(plasma2)], na.rm = TRUE)
  
  # Divide each protein by sum of that sample 
  ab <- function(x){x/column_sums_B}
  B3 <- apply(brain2[smpl_start:ncol(brain2)],1,ab)
  B3 <- t(B3)
  #add back the columns lost in calculations
  info_column_end <- (as.numeric(unlist(smpl_start)) - 1)
  B3 <- cbind(brain2[1:info_column_end], B3)
  
  # Divide each protein by sum of that sample 
  ac <- function(x){x/column_sums_C}
  C3 <- apply(CSF2[smpl_start:ncol(CSF2)],1,ac)
  C3 <- t(C3)
  #add back the columns lost in calculations
  info_column_end <- (as.numeric(unlist(smpl_start)) - 1)
  C3 <- cbind(CSF2[1:info_column_end], C3)
  
  # Divide each protein by sum of that sample 
  ap <- function(x){x/column_sums_P}
  P3 <- apply(plasma2[smpl_start:ncol(plasma2)],1,ap)
  P3 <- t(P3)
  #add back the columns lost in calculations
  info_column_end <- (as.numeric(unlist(smpl_start)) - 1)
  P3 <- cbind(plasma2[1:info_column_end], P3)
  
  # remove high abundant depleted proteins and keratin
  P3 <- P3 %>% filter(!Accession %in% high_list)
  B3 <- B3 %>% filter(!str_detect(Description, "KRT"))
  C3 <- C3 %>% filter(!str_detect(Description, "KRT"))
  P_KRT <- subset(P3, str_detect(Description, "KRT"))
  P3 <- P3 %>% filter(!str_detect(Description, "KRT"))
  
  
  #keep only if x% complete
  b90 <- ceiling(0.15*19) #brain 85% complete
  c90 <- ceiling(0.1*8) #CSF 90% complete
  p90 <- ceiling(0.3*20) #plasma 70% complete
  
  B4 <-  subset(B3, rowSums(is.na(B3[,c(8:17, 28:37)])) < b90 & rowSums(is.na(B3[,c(18:27, 38:46)])) < b90)
  B4 <-  subset(B4, B4$`Protein FDR Confidence: Combined` == "High" |  B4$`Protein FDR Confidence: Combined` == "Medium" )
  
  C4 <-  subset(C3, rowSums(is.na(C3[,8:15])) < c90 & rowSums(is.na(C3[,c(16:23)])) < c90)
  C4 <-  subset(C4, C4$`Protein FDR Confidence: Combined` == "High" |  C4$`Protein FDR Confidence: Combined` == "Medium" )
  
  P4 <-  subset(P3, rowSums(is.na(P3[,8:27])) < p90 & rowSums(is.na(P3[,28:47])) < p90)
  P4 <-  subset(P4, P4$`Protein FDR Confidence: Combined` == "High" |  P4$`Protein FDR Confidence: Combined` == "Medium" )
  
  #Subsetting groups
  #brain
  B_AD <- B4 %>% dplyr::select(matches(c("AD")))
  B_C <- B4 %>% dplyr::select(matches(c("C_[1-9]")))
  #CSF
  CSF_AD <- C4 %>% dplyr::select(matches(c("AD")))
  CSF_C <- C4 %>% dplyr::select(matches(c("C_V")))
  #plasma
  P_AD <- P4 %>% dplyr::select(matches(c("AD")))
  P_C <- P4 %>% dplyr::select(matches(c("C_[1-9]")))
  
  # for separation by gender
  Brain_Male_AD <- B4 %>% dplyr::select(matches(c("M_AD")))
  Brain_Male_Control <- B4 %>% dplyr::select(matches(c("M_C")))
  Brain_Female_AD <- B4 %>% dplyr::select(matches(c("F_AD")))
  Brain_Female_Control <- B4 %>% dplyr::select(matches(c("F_C")))
  
  CSF_Male_AD <- C4 %>% dplyr::select(matches(c("M_AD")))
  CSF_Male_Control <- C4 %>% dplyr::select(matches(c("M_C")))
  CSF_Female_AD <- C4 %>% dplyr::select(matches(c("F_AD")))
  CSF_Female_Control <- C4 %>% dplyr::select(matches(c("F_C")))
  
  Plasma_Male_AD <- P4 %>% dplyr::select(matches(c("M_AD")))
  Plasma_Male_Control <- P4 %>% dplyr::select(matches(c("M_C")))
  Plasma_Female_AD <- P4 %>% dplyr::select(matches(c("F_AD")))
  Plasma_Female_Control <- P4 %>% dplyr::select(matches(c("F_C")))
  
  # Cohort Average
  B_stats_df <- data.frame(B4$Accession, B4$`Protein FDR Confidence: Combined`, B4$Description)
  colnames(B_stats_df) <- (c("Accession", "Confidence", "Description"))
  B_stats_df$AD_avg <- rowMeans(B_AD, na.rm = TRUE)
  B_stats_df$Control_avg <- rowMeans(B_C, na.rm = TRUE)
  B_stats_df$M_AD_avg <- rowMeans(Brain_Male_AD, na.rm = TRUE)
  B_stats_df$M_Control_avg <- rowMeans(Brain_Male_Control, na.rm = TRUE)
  B_stats_df$F_AD_avg <- rowMeans(Brain_Female_AD, na.rm = TRUE)
  B_stats_df$F_Control_avg <- rowMeans(Brain_Female_Control, na.rm = TRUE)
  
  CSF_stats_df <- data.frame(C4$Accession, C4$`Protein FDR Confidence: Combined`, C4$Description)
  colnames(CSF_stats_df) <- (c("Accession", "Confidence", "Description"))
  CSF_stats_df$AD_avg <- rowMeans(CSF_AD, na.rm = TRUE)
  CSF_stats_df$Control_avg <- rowMeans(CSF_C, na.rm = TRUE)
  CSF_stats_df$M_AD_avg <- rowMeans(CSF_Male_AD, na.rm = TRUE)
  CSF_stats_df$M_Control_avg <- rowMeans(CSF_Male_Control, na.rm = TRUE)
  CSF_stats_df$F_AD_avg <- rowMeans(CSF_Female_AD, na.rm = TRUE)
  CSF_stats_df$F_Control_avg <- rowMeans(CSF_Female_Control, na.rm = TRUE)
  
  P_stats_df <- data.frame(P4$Accession, P4$`Protein FDR Confidence: Combined`, P4$Description)
  colnames(P_stats_df) <- (c("Accession", "Confidence", "Description"))
  P_stats_df$AD_avg <- rowMeans(P_AD, na.rm = TRUE)
  P_stats_df$Control_avg <- rowMeans(P_C, na.rm = TRUE)
  P_stats_df$M_AD_avg <- rowMeans(Plasma_Male_AD, na.rm = TRUE)
  P_stats_df$M_Control_avg <- rowMeans(Plasma_Male_Control, na.rm = TRUE)
  P_stats_df$F_AD_avg <- rowMeans(Plasma_Female_AD, na.rm = TRUE)
  P_stats_df$F_Control_avg <- rowMeans(Plasma_Female_Control, na.rm = TRUE)
  
  # Fold Change
  # Brain
  B_stats_df$AD_vs_C_fold_change <- (B_stats_df$AD_avg)/(B_stats_df$Control_avg)
  B_stats_df$M_AD_vs_C_fold_change <- (B_stats_df$M_AD_avg)/(B_stats_df$M_Control_avg)
  B_stats_df$F_AD_vs_C_fold_change <- (B_stats_df$F_AD_avg)/(B_stats_df$F_Control_avg)
  
  
  #CSF
  CSF_stats_df$AD_vs_C_fold_change <- (CSF_stats_df$AD_avg)/(CSF_stats_df$Control_avg)
  CSF_stats_df$M_AD_vs_C_fold_change <- (CSF_stats_df$M_AD_avg)/(CSF_stats_df$M_Control_avg)
  CSF_stats_df$F_AD_vs_C_fold_change <- (CSF_stats_df$F_AD_avg)/(CSF_stats_df$F_Control_avg)
  
  
  #Plasma
  P_stats_df$AD_vs_C_fold_change <- (P_stats_df$AD_avg)/(P_stats_df$Control_avg)
  P_stats_df$M_AD_vs_C_fold_change <- (P_stats_df$M_AD_avg)/(P_stats_df$M_Control_avg)
  P_stats_df$F_AD_vs_C_fold_change <- (P_stats_df$F_AD_avg)/(P_stats_df$F_Control_avg)
  
  
  
  #Log2 normalization
  #Brain
  B_stats_df$log2fc <- log(B_stats_df$AD_vs_C_fold_change, 2)
  B_stats_df$M_log2fc <- log(B_stats_df$M_AD_vs_C_fold_change, 2)
  B_stats_df$F_log2fc <- log(B_stats_df$F_AD_vs_C_fold_change, 2)
  
  #CSF
  CSF_stats_df$log2fc <- log(CSF_stats_df$AD_vs_C_fold_change, 2)
  CSF_stats_df$M_log2fc <- log(CSF_stats_df$M_AD_vs_C_fold_change, 2)
  CSF_stats_df$F_log2fc <- log(CSF_stats_df$F_AD_vs_C_fold_change, 2)
  
  #Plasma
  P_stats_df$log2fc <- log(P_stats_df$AD_vs_C_fold_change, 2)
  P_stats_df$M_log2fc <- log(P_stats_df$M_AD_vs_C_fold_change, 2)
  P_stats_df$F_log2fc <- log(P_stats_df$F_AD_vs_C_fold_change, 2)
  
  # P val calculations (Welch's t-test)
  # brain AD vs C
  B_stats_df$AD_C_pval <- sapply(1:nrow(B4), function(i) t.test(as.numeric(as.character(unlist(B_AD[i,1:20]))), as.numeric(as.character(unlist(B_C[i,1:19]))))[c("p.value")])
  B_stats_df$M_AD_C_pval <- sapply(1:nrow(B4), function(i) t.test(as.numeric(as.character(unlist(Brain_Male_AD[i,1:10]))), as.numeric(as.character(unlist(Brain_Male_Control[i,1:9]))))[c("p.value")])
  B_stats_df$F_AD_C_pval <- sapply(1:nrow(B4), function(i) t.test(as.numeric(as.character(unlist(Brain_Female_AD[i,1:10]))), as.numeric(as.character(unlist(Brain_Female_Control[i,1:9]))))[c("p.value")])
  
  
  #CSF
  CSF_stats_df$AD_C_pval <- sapply(1:nrow(C4), function(i) t.test(as.numeric(as.character(unlist(CSF_AD[i,1:8]))), as.numeric(as.character(unlist(CSF_C[i,1:8]))))[c("p.value")])
  CSF_stats_df$M_AD_C_pval <- sapply(1:nrow(C4), function(i) t.test(as.numeric(as.character(unlist(CSF_Male_AD[i,1:4]))), as.numeric(as.character(unlist(CSF_Male_Control[i,1:4]))))[c("p.value")])
  CSF_stats_df$F_AD_C_pval <- sapply(1:nrow(C4), function(i) t.test(as.numeric(as.character(unlist(CSF_Female_AD[i,1:4]))), as.numeric(as.character(unlist(CSF_Female_Control[i,1:4]))))[c("p.value")])
  
  #plasma
  P_stats_df$AD_C_pval <- sapply(1:nrow(P4), function(i) t.test(as.numeric(as.character(unlist(P_AD[i,1:20]))), as.numeric(as.character(unlist(P_C[i,1:20]))))[c("p.value")])
  P_stats_df$M_AD_C_pval <- sapply(1:nrow(P4), function(i) t.test(as.numeric(as.character(unlist(Plasma_Male_AD[i,1:10]))), as.numeric(as.character(unlist(Plasma_Male_Control[i,1:10]))))[c("p.value")])
  P_stats_df$F_AD_C_pval <- sapply(1:nrow(P4), function(i) t.test(as.numeric(as.character(unlist(Plasma_Female_AD[i,1:10]))), as.numeric(as.character(unlist(Plasma_Female_Control[i,1:10]))))[c("p.value")])
  
  
  # Select statistically significant for correction
  #Brain
  sig_B_AD <- subset(B_stats_df, B_stats_df$AD_C_pval < 0.05)
  sig_B_M_AD <- subset(B_stats_df, B_stats_df$M_AD_C_pval < 0.05)
  sig_B_F_AD <- subset(B_stats_df, B_stats_df$F_AD_C_pval < 0.05)
  
  #CSF
  sig_CSF_AD <- subset(CSF_stats_df, CSF_stats_df$AD_C_pval < 0.05)
  sig_CSF_M_AD <- subset(CSF_stats_df, CSF_stats_df$M_AD_C_pval < 0.05)
  sig_CSF_F_AD <- subset(CSF_stats_df, CSF_stats_df$F_AD_C_pval < 0.05)
  
  #Plasma
  sig_P_AD <- subset(P_stats_df, P_stats_df$AD_C_pval < 0.05)
  sig_P_M_AD <- subset(P_stats_df, P_stats_df$M_AD_C_pval < 0.05)
  sig_P_F_AD <- subset(P_stats_df, P_stats_df$F_AD_C_pval < 0.05) 
  
  #Venn Diagrams
  #Brain 
  a <- unlist(sig_B_AD[,1])
  b <- unlist(sig_B_M_AD[,1])
  c <- unlist(sig_B_F_AD[,1])
  
  venn <- list(`AD vs C` = a, 
               `Male` = b, 
               `Female` = c)
  
  ggvenn(venn, show_percentage = F, text_size = 7, set_name_size = 10, show_outside = "auto", fill_color = c("yellow3","steelblue","tomato3"), fill_alpha = 0.6, stroke_size = 1.5)
  
  #gold2, tomato, skyblue
  ggsave("plots/AD_bain_sig_venn_proteins_V2.svg", width = 7.5, height = 7.5, units = "in")
  
  v.table <- venn(venn)
  sec <- attributes(v.table)
  
  #extract sections
  brain_venn_all <- sec$intersections$`AD vs C:Male:Female`
  brain_venn_male <- sec$intersection$`Male`
  brain_venn_female <- sec$intersection$`Female`
  brain_venn_ad <- sec$intersection$`AD vs C`
  brain_venn_AD_and_male <- sec$intersection$`AD vs C:Male`
  brain_venn_AD_and_female <- sec$intersection$`AD vs C:Female`
  brain_venn_male_and_female <- sec$intersections$`Male:Female`
  
  write.xlsx(brain_venn_all,"processed_data/brain_sig_venn_sections.xlsx", sheetName = "all")
  write.xlsx(brain_venn_male, "processed_data/brain_sig_venn_sections.xlsx", sheetName = "male", append = TRUE)
  write.xlsx(brain_venn_female, "processed_data/brain_sig_venn_sections.xlsx", sheet = "female", append = TRUE)
  write.xlsx(brain_venn_ad, "processed_data/brain_sig_venn_sections.xlsx", sheet = "ad", append = TRUE)
  write.xlsx(brain_venn_AD_and_male, "processed_data/brain_sig_venn_sections.xlsx", sheet = "ad_and male", append = TRUE)
  write.xlsx(brain_venn_AD_and_female, "processed_data/brain_sig_venn_sections.xlsx", sheet = "ad_and female", append = TRUE)
  write.xlsx(brain_venn_male_and_female, "processed_data/brain_sig_venn_sections.xlsx", sheet = "male and female", append = TRUE)
  
  #CSF
  a <- unlist(sig_CSF_AD[,1])
  b <- unlist(sig_CSF_M_AD[,1])
  c <- unlist(sig_CSF_F_AD[,1])
  
  venn <- list(`AD vs C` = a, 
               `Male` = b, 
               `Female` = c)
  
  ggvenn(venn, show_percentage = F, text_size = 7, set_name_size = 10,
         show_outside = "auto", fill_color = c("yellow3","steelblue","tomato3"),
         fill_alpha = 0.6, stroke_size = 1.5)
  
  #gold2, tomato, skyblue
  ggsave("plots/AD_CSF_sig_venn_proteins_v2.svg", width = 7.5, height = 7.5, units = "in")
  
  v.table <- venn(venn)
  sec <- attributes(v.table)
  
  #extract sections
  csf_venn_all <- sec$intersections$`AD vs C:Male:Female`
  csf_venn_male <- sec$intersection$`Male`
  csf_venn_female <- sec$intersection$`Female`
  csf_venn_ad <- sec$intersection$`AD vs C`
  csf_venn_AD_and_male <- sec$intersection$`AD vs C:Male`
  csf_venn_AD_and_female <- sec$intersection$`AD vs C:Female`
  csf_venn_male_and_female <- sec$intersections$`Male:Female`
  
  write.xlsx(csf_venn_all,"processed_data/csf_sig_venn_sections.xlsx", sheetName = "all")
  write.xlsx(csf_venn_male, "processed_data/csf_sig_venn_sections.xlsx", sheetName = "male", append = TRUE)
  write.xlsx(csf_venn_female, "processed_data/csf_sig_venn_sections.xlsx", sheet = "female", append = TRUE)
  write.xlsx(csf_venn_ad, "processed_data/csf_sig_venn_sections.xlsx", sheet = "ad", append = TRUE)
  write.xlsx(csf_venn_AD_and_male, "processed_data/csf_sig_venn_sections.xlsx", sheet = "ad_and male", append = TRUE)
  write.xlsx(csf_venn_AD_and_female, "processed_data/csf_sig_venn_sections.xlsx", sheet = "ad_and female", append = TRUE)
  write.xlsx(csf_venn_male_and_female, "processed_data/csf_sig_venn_sections.xlsx", sheet = "male and female", append = TRUE)
  
  #Plasma
  a <- unlist(sig_P_AD[,1])
  b <- unlist(sig_P_M_AD[,1])
  c <- unlist(sig_P_F_AD[,1])
  
  venn <- list(`AD vs C` = a, 
               `Male` = b, 
               `Female` = c)
  
  ggvenn(venn, show_percentage = F, text_size = 7, set_name_size = 10, show_outside = "auto", fill_color = c("yellow3","steelblue","tomato3"), fill_alpha = 0.6, stroke_size = 1.5)
  
  #gold2, tomato, skyblue
  ggsave("plots/AD_plasma_sig_venn_proteins.svg", width = 7.5, height = 7.5, units = "in")
  
  v.table <- venn(venn)
  sec <- attributes(v.table)
  
  #extract sections
  plasma_venn_all <- sec$intersections$`AD vs C:Male:Female`
  plasma_venn_male <- sec$intersection$`Male`
  plasma_venn_female <- sec$intersection$`Female`
  plasma_venn_ad <- sec$intersection$`AD vs C`
  plasma_venn_AD_and_male <- sec$intersection$`AD vs C:Male`
  plasma_venn_AD_and_female <- sec$intersection$`AD vs C:Female`
  plasma_venn_male_and_female <- sec$intersections$`Male:Female`
  
  write.xlsx(plasma_venn_all,"processed_data/plasma_sig_venn_sections.xlsx", sheetName = "all")
  write.xlsx(plasma_venn_male, "processed_data/plasma_sig_venn_sections.xlsx", sheetName = "male", append = TRUE)
  write.xlsx(plasma_venn_female, "processed_data/plasma_sig_venn_sections.xlsx", sheet = "female", append = TRUE)
  write.xlsx(plasma_venn_ad, "processed_data/plasma_sig_venn_sections.xlsx", sheet = "ad", append = TRUE)
  write.xlsx(plasma_venn_AD_and_male, "processed_data/plasma_sig_venn_sections.xlsx", sheet = "ad_and male", append = TRUE)
  write.xlsx(plasma_venn_AD_and_female, "processed_data/plasma_sig_venn_sections.xlsx", sheet = "ad_and female", append = TRUE)
  write.xlsx(plasma_venn_male_and_female, "processed_data/plasma_sig_venn_sections.xlsx", sheet = "male and female", append = TRUE)
  
  #compare statistically significant across sample type
  
  #3 sample types
  a <- unlist(sig_B_AD[,1])
  b <- unlist(sig_CSF_AD[,1])
  c <- unlist(sig_P_AD[,1])
  
  venn <- list(`Brain` = a, 
               `CSF` = b, 
               `Plasma` = c)
  
  ggvenn(venn, show_percentage = F, text_size = 7, set_name_size = 10, show_outside = "auto", fill_color = c("yellow3","steelblue","tomato3"), fill_alpha = 0.6, stroke_size = 1.5)
  
  #gold2, tomato, skyblue
  ggsave("plots/AD_3_sample_sig_venn_proteins.svg", width = 7.5, height = 7.5, units = "in")
  
  v.table <- venn(venn)
  sec <- attributes(v.table)
  
  #extract sections
  sig3_venn_all <- sec$intersections$`Brain:CSF:Plasma`
  sig3_venn_brain <- sec$intersection$Brain
  sig3_venn_CSF <- sec$intersection$CSF
  sig3_venn_plasma <- sec$intersection$Plasma
  sig3_venn_CSF_plasma <- sec$intersection$`CSF:Plasma`
  sig3_venn_CSF_brain <- sec$intersection$`Brain:CSF`
  sig3_venn_brain_plasma <- sec$intersections$`Brain:Plasma`
  
  write.xlsx(sig3_venn_all,"processed_data/sig3_venn_sections.xlsx", sheetName = "all")
  write.xlsx(sig3_venn_brain, "processed_data/sig3_venn_sections.xlsx", sheetName = "brain", append = TRUE)
  write.xlsx(sig3_venn_CSF, "processed_data/sig3_venn_sections.xlsx", sheet = "CSF", append = TRUE)
  write.xlsx(sig3_venn_plasma, "processed_data/sig3_venn_sections.xlsx", sheet = "plasma", append = TRUE)
  write.xlsx(sig3_venn_CSF_plasma, "processed_data/sig3_venn_sections.xlsx", sheet = "CSF_plasma", append = TRUE)
  write.xlsx(sig3_venn_CSF_brain, "processed_data/sig3_venn_sections.xlsx", sheet = "CSF_brain", append = TRUE)
  write.xlsx(sig3_venn_brain_plasma, "processed_data/sig3_venn_sections.xlsx", sheet = "brain_plasma", append = TRUE)
  
  
  B5 <- B4[,c(1:3, 8:46)]
  C5 <- C4[, c(1:3, 8:23)]
  P5 <- P4[, c(1:3, 8:47)]
  
  B_age_sex <- read_excel("data/age_gender_ID.xlsx", sheet = "brain")
  C_age_sex <- read_excel("data/age_gender_ID.xlsx", sheet = "csf")
  P_age_sex <- read_excel("data/age_gender_ID.xlsx", sheet = "plasma")
  
  long_brain <- B5 %>% pivot_longer(cols = c(F_AD_1:M_C_10), 
                                    names_to = "Sample",
                                    values_to = "Abundance")
  long_CSF <- C5 %>% pivot_longer(cols = c(M_AD_327:M_C_V85500), 
                                  names_to = "Sample", 
                                  values_to = "Abundance")
  long_plasma <- P5 %>% pivot_longer(cols = c(F_C_25176:F_AD_84888), 
                                     names_to = "Sample", 
                                     values_to = "Abundance")
  long_brain <- merge(long_brain, B_age_sex, by = "Sample")
  long_CSF <- merge(long_CSF, C_age_sex, by = "Sample")
  long_plasma <- merge(long_plasma, P_age_sex, by = "Sample")
  
  #save long data as CSV
  write.csv(long_brain, "processed_data/long_brain_Sept_2024.csv")
  write.csv(long_CSF, "processed_data/long_CSF_Sept_2024.csv")
  write.csv(long_plasma, "processed_data/long_plasma_Sept_2024.csv")
  
  ####### RESULTS #######
  write.xlsx(B_stats_df,"processed_data/AD_proteomics_R_results_Brain(Sept_2024).xlsx", sheetName = "stats")
  write.xlsx(sig_B_AD,"processed_data/AD_proteomics_R_results_Brain(Sept_2024).xlsx", sheetName = "sig", append = TRUE)
  write.xlsx(sig_B_M_AD,"processed_data/AD_proteomics_R_results_Brain(Sept_2024).xlsx", sheetName = "Male", append = TRUE)
  write.xlsx(sig_B_F_AD,"processed_data/AD_proteomics_R_results_Brain(Sept_2024).xlsx", sheetName = "Female", append = TRUE)
  
  write.xlsx(CSF_stats_df,"processed_data/AD_proteomics_R_results_CSF(Sept_2024).xlsx", sheetName = "stats")
  write.xlsx(sig_CSF_AD,"processed_data/AD_proteomics_R_results_CSF(Sept_2024).xlsx", sheetName = "sig", append = TRUE)
  write.xlsx(sig_CSF_M_AD,"processed_data/AD_proteomics_R_results_CSF(Sept_2024).xlsx", sheetName = "Male", append = TRUE)
  write.xlsx(sig_CSF_F_AD,"processed_data/AD_proteomics_R_results_CSF(Sept_2024).xlsx", sheetName = "Female", append = TRUE)
  
  write.xlsx(P_stats_df,"processed_data/AD_proteomics_R_results_Plasma(Sept_2024).xlsx", sheetName = "stats")
  write.xlsx(sig_P_AD,"processed_data/AD_proteomics_R_results_Plasma(Sept_2024).xlsx", sheetName = "sig", append = TRUE)
  write.xlsx(sig_P_M_AD,"processed_data/AD_proteomics_R_results_Plasma(Sept_2024).xlsx", sheetName = "Male", append = TRUE)
  write.xlsx(sig_P_F_AD,"processed_data/AD_proteomics_R_results_Plasma(Sept_2024).xlsx", sheetName = "Female", append = TRUE)
}

#### Principal Component Analysis ####
# Plasma 
# AD vs C

###################################

# 2D PCA plot
#packages
long_brain <- read.csv("processed_data/long_brain_Sept_2024.csv")
long_CSF <- read.csv("processed_data/long_CSF_Sept_2024.csv")
long_plasma <- read.csv("processed_data/long_plasma_Sept_2024.csv")

long_brain <- long_brain %>% dplyr::select(c(2,4,6))
long_CSF <- long_CSF %>% dplyr::select(c(2,4,6))
long_plasma <- long_plasma %>% dplyr::select(c(2,4,6))


brain_wide <- long_brain %>% 
  pivot_wider(names_from = Sample, 
              values_from = Abundance)
csf_wide <- long_CSF %>% 
  pivot_wider(names_from = Sample, 
              values_from = Abundance)
plasma_wide <- long_plasma %>% 
  pivot_wider(names_from = Sample, 
              values_from = Abundance)



B_no_na <- na.omit(brain_wide)
C_no_na <- na.omit(csf_wide)
P_no_na <- na.omit(plasma_wide)

norm_pca <-C_no_na[,c(2:17)] 

norm_pca <- data.frame(t(norm_pca))
norm_pca$Sample <- c(rep("C (N=8)", 8), rep("AD (N=8)", 8))
norm_pca <- norm_pca %>% relocate(Sample)
data.pca <- prcomp(norm_pca[,c(2:ncol(norm_pca))], 
                   center = TRUE, 
                   scale. = TRUE)

pca_ad <- fviz_pca_ind(data.pca,
                       geom.ind = c("point"), # show points only (but not "text")
                       pointsize = 2,
                       col.ind = norm_pca$Sample, # color by groups
                       palette = c("purple3", "gray30"),
                       legend.title = "Groups", 
                       addEllipses = T,
                       ellipse.level = 0.95,
                       title = "PCA (CSF)", 
                       mean.point = F)

pca_ad <- pca_ad +theme_light()+ theme(legend.position = "bottom")+theme(plot.title=element_text(hjust=0.5))
pca_ad
ggsave("plots/new_pca_AD_vs_C_proteins_CS.svg", width = 5, height = 5, units = "in") 

###
# Plasma 
# Female (AD vs C)
pca_female <- C_no_na %>% dplyr::select(matches(c("F_C","F_AD")))
norm_pca <- data.frame(t(pca_female))
norm_pca$Sample <- c(rep("FC (N=4)", 4), rep("FAD (N=4)", 4))
norm_pca <- norm_pca %>% relocate(Sample)
data.pca <- prcomp(norm_pca[,c(2:ncol(norm_pca))], 
                   center = TRUE, 
                   scale. = TRUE)

pca_female <- fviz_pca_ind(data.pca,
                           geom.ind = c("point"), # show points only (but not "text")
                           pointsize = 2,
                           col.ind = norm_pca$Sample, # color by groups
                           palette = c("magenta", "pink3"),
                           legend.title = "Groups", 
                           addEllipses = T,
                           ellipse.level = 0.95,
                           title = "PCA of CSF Proteins (Female)", 
                           mean.point = F)
pca_female <- pca_female +theme_light()+ theme(legend.position = "bottom")+theme(plot.title=element_text(hjust=0.5))
pca_female
ggsave("plots/new_pca_AD_female_proteins_CSF.svg", width = 5, height = 5, units = "in") 

###
# Plasma 
# Male (AD vs C)
pca_male <-  C_no_na %>% dplyr::select(matches(c("M_C","M_AD")))
norm_pca <- data.frame(t(pca_male))
norm_pca$Sample <- c(rep("MC (N=4)", 4), rep("MAD (N=4)", 4))
norm_pca <- norm_pca %>% relocate(Sample)
data.pca <- prcomp(norm_pca[,c(2:ncol(norm_pca))], 
                   center = TRUE, 
                   scale. = TRUE)

pca_male <- fviz_pca_ind(data.pca,
                         geom.ind = c("point"), # show points only (but not "text")
                         pointsize = 2,
                         col.ind = norm_pca$Sample, # color by groups
                         palette = c("blue3", "steelblue"),
                         legend.title = "Groups", 
                         addEllipses = T,
                         ellipse.level = 0.95,
                         title = "PCA of CSF Proteins (Male)", 
                         mean.point = F)

pca_male <- pca_male +theme_light()+ theme(legend.position = "bottom")+theme(plot.title=element_text(hjust=0.5))
pca_male
ggsave("plots/pca_AD_male_proteins_male_CSF.svg", width = 5, height = 5, units = "in") 

###################################

#### Volcano plots ####

#select data
clean_data <- B_stats_df

#User Input Customization
graph_title <-  as.character(readline(prompt = "volcano title? :" ))

{
  clean_data$diffexpressed <- "NO"
  clean_data$diffexpressed[clean_data$log2fc > 1 & clean_data$AD_C_pval < 0.05] <- "UP"
  clean_data$diffexpressed[clean_data$log2fc < -1 & clean_data$AD_C_pval < 0.05] <- "DOWN"
  
  up_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "UP")))
  down_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "DOWN")))
  ns <- as.numeric(sum(str_count(clean_data$diffexpressed, "NO")))
}

#Total proteins
total_prots <- (ns+down_reg+up_reg)

ymax <- -log(min(unlist(clean_data$AD_C_pval)), 10)
xmin <- min(clean_data$log2fc)
xmax <- max(clean_data$log2fc)

annotation <- data.frame(
  x = c((xmin*0.75),(xmax*0.75), -Inf),
  y = c((ymax*0.5), (ymax*0.5), Inf),
  label = c(down_reg, up_reg, paste0("Total Proteins: ", total_prots)))

#Theme Set
theme_set(theme_bw())

v <- ggplot(data = clean_data, aes(x = log2fc, y = -log10(unlist(AD_C_pval)), col = diffexpressed, )) + 
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = c(1.3), col = "black", linetype = "dashed", linewidth = 1) +
  geom_point(size = 1.5) + 
  #geom_text_repel(label = B_stats_df$Accession, max.overlaps = Inf) +
  
  geom_label(data=annotation, aes(x=x, y=y, label=label, hjust="inward", vjust="inward"),
             color="black", 
             size=4, angle=0, fontface="bold") + 
  ggtitle(graph_title) + 
  labs(col = "Glycans") + labs(x = expression(log[2]*"fold change"), y = expression(-log[10]*"p-value")) +
  
  scale_color_manual(values = c("blue4", "grey60", "orange3"), 
                     labels = c("Downregulated", "Outside Statistical Parameters", "Upregulated"))

### optional colors for volcano:
### down(purple3, red3, blue4, firebrick)
### up (yellow3, green3, orange3, forestgreen)

v <- v+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 12))+theme(legend.text = element_text(size = 10))+guides(color = guide_legend(override.aes = list(size = 4)))+theme(legend.title = element_blank())+theme(legend.position = "bottom")+theme(plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

v

ggsave("plots/volcano_AD_brain.svg", width = 7, height = 6, units = "in")

###############################

#select data
clean_data <- CSF_stats_df

#User Input Customization
graph_title <-  as.character(readline(prompt = "volcano title? :" ))

{
  clean_data$diffexpressed <- "NO"
  clean_data$diffexpressed[clean_data$log2fc > 1 & clean_data$AD_C_pval < 0.05] <- "UP"
  clean_data$diffexpressed[clean_data$log2fc < -1 & clean_data$AD_C_pval < 0.05] <- "DOWN"
  
  up_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "UP")))
  down_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "DOWN")))
  ns <- as.numeric(sum(str_count(clean_data$diffexpressed, "NO")))
}

#Total proteins
total_prots <- (ns+down_reg+up_reg)

ymax <- -log(min(unlist(clean_data$AD_C_pval)), 10)
xmin <- min(clean_data$log2fc)
xmax <- max(clean_data$log2fc)

annotation <- data.frame(
  x = c((xmin*0.75),(xmax*0.75), -Inf),
  y = c((ymax*0.5), (ymax*0.5), Inf),
  label = c(down_reg, up_reg, paste0("Total Proteins: ", total_prots)))

#Theme Set
theme_set(theme_bw())

v <- ggplot(data = clean_data, aes(x = log2fc, y = -log10(unlist(AD_C_pval)), col = diffexpressed, )) + 
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = c(1.3), col = "black", linetype = "dashed", linewidth = 1) +
  geom_point(size = 1.5) + 
  #geom_text_repel(label = B_stats_df$Accession, max.overlaps = Inf) +
  
  geom_label(data=annotation, aes(x=x, y=y, label=label, hjust="inward", vjust="inward"),
             color="black", 
             size=4, angle=0, fontface="bold") + 
  ggtitle(graph_title) + 
  labs(col = "Glycans") + labs(x = expression(log[2]*"fold change"), y = expression(-log[10]*"p-value")) +
  
  scale_color_manual(values = c("blue4", "grey60", "orange3"), 
                     labels = c("Downregulated", "Outside Statistical Parameters", "Upregulated"))

### optional colors for volcano:
### down(purple3, red3, blue4, firebrick)
### up (yellow3, green3, orange3, forestgreen)

v <- v+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 12))+theme(legend.text = element_text(size = 10))+guides(color = guide_legend(override.aes = list(size = 4)))+theme(legend.title = element_blank())+theme(legend.position = "bottom")+theme(plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

v

ggsave("plots/volcano_AD_CSF.svg", width = 7, height = 6, units = "in")

#######################

#select data
clean_data <- P_stats_df

#User Input Customization
graph_title <-  as.character(readline(prompt = "volcano title? :" ))

{
  clean_data$diffexpressed <- "NO"
  clean_data$diffexpressed[clean_data$log2fc > 1 & clean_data$AD_C_pval < 0.05] <- "UP"
  clean_data$diffexpressed[clean_data$log2fc < -1 & clean_data$AD_C_pval < 0.05] <- "DOWN"
  
  up_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "UP")))
  down_reg <- as.numeric(sum(str_count(clean_data$diffexpressed, "DOWN")))
  ns <- as.numeric(sum(str_count(clean_data$diffexpressed, "NO")))
}

#Total proteins
total_prots <- (ns+down_reg+up_reg)

ymax <- -log(min(unlist(clean_data$AD_C_pval)), 10)
xmin <- min(clean_data$log2fc)
xmax <- max(clean_data$log2fc)

annotation <- data.frame(
  x = c((xmin*0.75),(xmax*0.75), -Inf),
  y = c((ymax*0.5), (ymax*0.5), Inf),
  label = c(down_reg, up_reg, paste0("Total Proteins: ", total_prots)))

#Theme Set
theme_set(theme_bw())

v <- ggplot(data = clean_data, aes(x = log2fc, y = -log10(unlist(AD_C_pval)), col = diffexpressed, )) + 
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed", linewidth = 1) +
  geom_hline(yintercept = c(1.3), col = "black", linetype = "dashed", linewidth = 1) +
  geom_point(size = 1.5) + 
  #geom_text_repel(label = B_stats_df$Accession, max.overlaps = Inf) +
  
  geom_label(data=annotation, aes(x=x, y=y, label=label, hjust="inward", vjust="inward"),
             color="black", 
             size=4, angle=0, fontface="bold") + 
  ggtitle(graph_title) + 
  labs(col = "Glycans") + labs(x = expression(log[2]*"fold change"), y = expression(-log[10]*"p-value")) +
  
  scale_color_manual(values = c("blue4", "grey60", "orange3"), 
                     labels = c("Downregulated", "Outside Statistical Parameters", "Upregulated"))

### optional colors for volcano:
### down(purple3, red3, blue4, firebrick)
### up (yellow3, green3, orange3, forestgreen)

v <- v+theme(axis.title = element_text(size = 14))+theme(axis.text = element_text(size = 12))+theme(legend.text = element_text(size = 10))+guides(color = guide_legend(override.aes = list(size = 4)))+theme(legend.title = element_blank())+theme(legend.position = "bottom")+theme(plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5))

v

ggsave("plots/volcano_AD_plasma.svg", width = 7, height = 6, units = "in")

#### Heat maps ####

#heatmap data set up

wide_brain <- pivot_wider(long_brain,
                          id_cols = "Accession",
                          names_from = "Sample", 
                          values_from = "Abundance")
wide_CSF <- pivot_wider(long_CSF,
                        id_cols = "Accession",
                        names_from = "Sample", 
                        values_from = "Abundance")
wide_plasma <- pivot_wider(long_plasma,
                           id_cols = "Accession",
                           names_from = "Sample", 
                           values_from = "Abundance")

#heat_data <- subset(long_brain, brain$Accession %in% long_brain$Accession)
wide_brain2 <- na.omit(wide_brain)
wide_CSF2 <- na.omit(wide_CSF)
wide_plasma2 <- na.omit(wide_plasma)

#subset only up or down + sig
UDS_brain <- subset(clean_data, clean_data$diffexpressed == "UP" | clean_data$diffexpressed =="DOWN")
UDS_brain_proteins <- UDS_brain$Accession

UDS_CSF <- subset(clean_data, clean_data$diffexpressed == "UP" | clean_data$diffexpressed =="DOWN")
UDS_CSF_proteins <- UDS_CSF$Accession

UDS_plasma <- subset(clean_data, clean_data$diffexpressed == "UP" | clean_data$diffexpressed =="DOWN")
UDS_plasma_proteins <- UDS_plasma$Accession

write.xlsx(UDS_brain_proteins, "processed_data/Up_down_sig AD proteins.xlsx", sheetName = "brain")
write.xlsx(UDS_CSF_proteins, "processed_data/Up_down_sig AD proteins.xlsx", sheetName = "csf", append = T)
write.xlsx(UDS_plasma_proteins, "processed_data/Up_down_sig AD proteins.xlsx", sheetName = "plasma", append = T)

heat_wide_brain <- subset(wide_brain, wide_brain$Accession%in%UDS_brain_proteins)
heat_wide_csf <- subset(wide_CSF, wide_CSF$Accession%in%UDS_CSF_proteins)
heat_wide_plasma <- subset(wide_plasma, wide_plasma$Accession%in%UDS_plasma_proteins)

#change for each set
heat_data <- heat_wide_plasma
#remove row names
proteins <- heat_data$Accession
proteins <- data.frame(proteins)
heat_data2 <- heat_data[,-1]
heat_map_data <- as.matrix(heat_data2)
row.names(heat_map_data) <- heat_data$Accession

#colr <- colorRampPalette(c("dodgerblue", "white", "orange"))(100)
colr <- colorRampPalette(c("dodgerblue", "white", "orange"))(100)
#dodgerblue
#gray0
#green3
#
# heatmap
h <- pheatmap(
  heat_map_data  ,
  color = colr,
  scale = "row",
  angle_col = "315",
  cluster_rows = T,
  cluster_cols = F,
  main = "Heatmap of Statistically Significant Up and Downregulated Proteins (Plasma)",
  fontsize = 16, 
  show_rownames = T,
  show_colnames = T,
  fontsize_row = 10,
  na_col = "gray50"
)
h
ggsave(here("plots/heatmap_AD_plasma_up_down_proteins.svg"), height = 6, width = 18, units = "in", plot = h)

#### Box plots ####
{
  #set_here("C:\\Users\\Andy Bennett\\Desktop\\R\\AD Proteomics\\manuscript_data\\")
  
  long_brain <- read.csv("long_brain.csv")
  
  AmPP <- subset(long_brain, long_brain$Accession == "P05067")
  
  b <- ggplot(AmPP, aes(x = AGE, y=  Abundance, color = Disease_State))+ 
    geom_point(alpha = 0.7, size = 2.5) + scale_color_manual(values = c("purple3","gray40"))                                  
  b  
  
  # add linear model
  b + geom_smooth(method = "lm")
  #?geom_smooth
  
  long_brain <- read.csv("processed_data/long_brain_Sept_2024.csv")
  long_csf <- read.csv("processed_data/long_CSF_Sept_2024.csv")
  long_plasma <- read.csv("processed_data/long_plasma_Sept_2024.csv")
  in_all <- read.xlsx("processed_data/sig3_venn_sections.xlsx", sheetName = "all")
  in_all <- in_all$x
  
  genes <- c("DKK3", "HSPG2", "IGHA2")
  # Combine the lists into a data frame
  IDs <- data.frame(accession = unlist(in_all), gene = unlist(genes))
  IDs
  
  long_plasma$Disease_State <- factor(long_plasma$Disease_State, levels = c("Control", "AD"))
  long_csf$Disease_State <- factor(long_csf$Disease_State, levels = c("Control", "AD"))
  long_brain$Disease_State <- factor(long_brain$Disease_State, levels = c("Control", "AD"))
  #long_plasma2 <- na.omit(long_plasma)
  # 3 found in all
  theme_set(theme_bw())
  for (protein in in_all) { 
    box1 <- subset(long_plasma, long_plasma$Accession == protein)
    box1$gene <- IDs$gene[IDs$accession == protein]
    # Plotting the BoxPlot
    {
      ymax <- max(box1$Abundance, na.rm = T) * 1.2
      ggplot(box1, aes(x = `Disease_State`, y = Abundance, fill = `Disease_State`)) +
        geom_boxplot(alpha=0.7, width = 0.6, lwd = 1, outlier.shape = 1, 
                     outlier.size = 1.6, outlier.alpha = 0.7, outlier.color='black', na.rm = TRUE) +
        xlab(NULL) +
        ylab("Relative Abundance") +
        scale_fill_manual(values = c("AD" = "purple3", "Control" = "gray40"))  +
        theme(axis.title.x = element_text(size = 20, face = "bold"),
              axis.title.y = element_text(size = 24, face = "bold"),
              axis.text.x = element_text(size = 28, face = "bold", color = "black"),
              axis.text.y = element_text(size = 24),
              panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA)) +
        
        theme(
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          #axis.line = element_line(color = "black", size = 1.1),
          plot.title = element_text(hjust = 0.5, size = 32, face = "bold")) +
        stat_compare_means(comparisons = list(c("AD", "Control")), method = "t.test", 
                           label = "p.signif",label.x = 1, size = 18, 
                           face = "bold", label.y = ymax*0.85, na.rm = TRUE)+ 
        scale_y_continuous(limits = c(0, ymax)) + labs(title = paste0(box1$gene," (", protein,")", " Plasma")) #+ 
      #theme(plot.title = element_markdown())
      
      
    }
    
    ggsave(here("processed_data/new_box",paste0("bold_AD_plasma_boxplot_",protein, ".svg")), height = 7, width = 7, units = "in", dpi = 600)
    
  }   
  
  ??element_markdown
  test <- subset(long_plasma2, long_plasma2$Accession == in_all)
  three <- unique(test[,4:5])
  
  table3 <- gt(three)
  table3
  gtsave(table3, "plots/3_FIA.docx")
  
  
  test2 <- subset(long_plasma, long_plasma$Accession == bp)
  fifteen <- unique(test2[,4:5])
  
  table15 <- gt(fifteen)
  table15
  gtsave(table15, "plots/15_bp.docx")
}

#### ROC analysis ####
{
  #brain <- read_excel("manuscript_data/NEW_AD_proteomics_R_results_Brain(April_2024).xlsx", sheet = "sig")
  #csf <- read_excel("manuscript_data/NEW_AD_proteomics_R_results_CSF(April_2024).xlsx", sheet = "sig")
  #plasma <- read_excel("manuscript_data/NEW_AD_proteomics_R_results_Plasma(April_2024).xlsx", sheet = "sig")
  
  #long_brain <- read.csv("manuscript_data/long_brain.csv")
  #long_csf <- read.csv("manuscript_data/long_csf.csv")
  long_plasma <- read.csv("processed_data/long_plasma_Sept_2024.csv")

  long_plasma2 <- long_plasma %>% 
    mutate(target = case_when(
      str_detect(Sample, "AD") ~ "1",
      TRUE ~ "0"
    ))
  
  unique_accessions <- long_plasma2 %>% 
    distinct(Accession)
  
  ## Just arranged them so 0s and 1s are together, doesnt cause any problems if skipped
  roc_data <- long_plasma2 %>% 
    filter(Accession == unique_accessions$Accession[1]) %>% 
    arrange(target)
  
  # Just a test
  roc_test <- pROC::roc(roc_data$target, roc_data$Abundance)
  
  ### ROC loop ###
  roc_df <- data.frame() # empty dataframe to store the roc results
  
  #select data for each ROC
  #3 in all
  in_all_3 <- read.xlsx("processed_data/sig3_venn_sections.xlsx", sheetName = "all")
  in_all_3 <- in_all_3[,2]
  
  #brain_plasma
  bp <- read.xlsx("processed_data/sig3_venn_sections.xlsx", sheetName = "brain_plasma")
  bp <- bp[,2]
  #male plasma unique
  male_prots <- read.xlsx("processed_data/plasma_sig_venn_sections.xlsx", sheetName = "male")
  male_prots <- male_prots[,2]
  
  #female plasma unique
  female_prots <- read.xlsx("processed_data/plasma_sig_venn_sections.xlsx", sheetName = "female")
  female_prots <- female_prots[,2]
  
  # Loop through each unique accession
  for (accession in bp) {
    roc_data <- long_plasma2 %>% 
      filter(Accession == accession) %>% 
      arrange(target)
    
    roc_test <- pROC::roc(roc_data$target, roc_data$Abundance) ## perform ROC
    sensitivity <- roc_test$sensitivities
    specificity <- 1 - roc_test$specificities
    auc <- rep(roc_test$auc, length(sensitivity)) # create a auc cplumn
    accessions <- rep(accession, length(sensitivity))
    
    roc_df <- rbind(roc_df, data.frame(accessions, sensitivity, specificity, auc))
    roc_df <- roc_df %>% 
      arrange(sensitivity, specificity, accessions)
  }

  ##############################################
  
  ## Binary Logistic Regression
  long_plasma3 <- long_plasma2 %>% dplyr::select("Accession", "Abundance", "target", "Sample")
  wide_plasma <- pivot_wider(long_plasma3,
                             names_from = Accession, 
                             values_from = Abundance)
  target <- wide_plasma$target
  wide_plasma <- wide_plasma %>% dplyr::select(-c("Sample"))
  
  wide_plasma2 <- wide_plasma %>% dplyr::select(target, contains(male_prots)) %>% 
    arrange(target)
  

  
  ##### Heading 1 #####
  wide_plasma2 <- wide_plasma2 %>% mutate_all(~as.numeric(as.character(.)))
  
  svm_model <- svm(target ~ ., data = wide_plasma2, probability = TRUE, na.action = na.exclude, kernel = "radial")

  png_d1_probabilities <- predict(svm_model, type = "response")
  
  print(png_d1_probabilities)
  target <- wide_plasma2$target
  png_d1_combined_roc <- roc(target, png_d1_probabilities)
  
  auc_value <- auc(png_d1_combined_roc)
  
  png_d1_combined_roc_data <- data.frame(
    
    accessions = rep("Combined", length(png_d1_combined_roc$sensitivities)),
    
    sensitivity = png_d1_combined_roc$sensitivities,
    specificity = 1 - png_d1_combined_roc$specificities, 
    auc = rep(auc_value, length(png_d1_combined_roc$sensitivities)),
    specificity2 = png_d1_combined_roc$specificities
    
  )
  # IMPORTANT!!!
  #r bind to roc_df
  ########################################
  roc_2 <- png_d1_combined_roc_data[,-5]
  #roc_2 <- roc_2 %>% 
  #arrange(sensitivity, specificity, accessions)
  
  ###compute average before Combined is added!!!
  avg_AUC <- signif(mean(roc_df$auc), digits = 2)
  roc_df <- rbind(roc_df, roc_2)
  
  
  ## Merge the accessions with the AUC values for the legend
  roc_df$auc <- signif(roc_df$auc, digits = 2)
  roc_df$legend <- paste(roc_df$accessions, " (AUC: ", roc_df$auc, ")", sep = "")
  
  roc_df_top <- roc_df %>% 
    filter(auc > 0.49) %>%
    arrange(sensitivity, specificity, accessions)## filter for AUC > 0.95
  
  r <- ggplot(roc_df_top, aes(x = specificity, y = sensitivity, color = legend)) +
    geom_step(size = 2) +
    geom_point(size = .05) +
    theme_pubr() +
    labs(title = "ROC Curve for Statistically Significant Proteins",
         subtitle = "Shared Between Brain Tissue and Plasma Samples",
         x = "1 - Specificity",
         y = "Sensitivity", 
         color = paste0("Accession ID (AUC)")) + scale_colour_manual(values=unname(alphabet()))+ 
    theme(axis.title.x = element_text(size = 18, face = "bold"),
          axis.title.y = element_text(size = 18, face = "bold"),
          axis.text.x = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 16, face = "bold")) +
    #panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA)) +
    #ggtitle("ROC Curve for Plasma Proteomics") +
    theme(
      legend.position = c(.75, .25),
      legend.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      #legend.text = element_text("Proteins with AUC"),
      #panel.grid.major = element_blank(),
      #panel.grid.minor = element_blank(),
      #panel.background = element_blank(),
      #axis.line = element_line(color = "black", size = 1.1),
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"), 
      plot.subtitle = element_text(hjust = 0.5, size = 16, face = "bold")) +
    geom_abline(intercept = 0, slope = 1, colour = "gray50", alpha = 0.5,linewidth = 1)
  #labs(color = "Protein (AUC)")

  # removing padding from roc curve
  r
  
  # moves axis to be exactly on plot
  r1 <- r+scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))
  r1
  
  #adjusts space around the whole figure
  r2 <- r1+ theme(plot.margin = unit(c(1,1,1,1), "cm"))
  r2
  r2 + guides(color = guide_legend(ncol = 2))
  
  ###############
  ggsave(here("plots/protein_roc_brain_and_plasma_final.svg"), width = 9, height = 9, units = "in")
  
  #Plasma
  a <- unlist(plasma[,2])
  b <- unlist(plasma_male[,2])
  c <- unlist(plasma_female[,2])
  
  venn <- list(`AD vs C` = a, 
               `Male` = b, 
               `Female` = c)
  
  ggvenn(venn, show_percentage = F, text_size = 7, set_name_size = 10, show_outside = "auto", fill_color = c("yellow3","steelblue","tomato3"), fill_alpha = 0.6, stroke_size = 1.5)
  
  #gold2, tomato, skyblue
  ggsave("plots/AD_sig_venn_proteins.svg", width = 7.5, height = 7.5, units = "in")
  
  v.table <- venn(venn)
  sec <- attributes(v.table)
  
  #extract sections
  plasma_in_all <- sec$intersections$`AD vs C:Male:Female`
  plasma_venn_male <- sec$intersections$Male
  plasma_venn_female <- sec$intersections$Female
}

#### GT tables ####
# 3 in all
gt3 <- subset(long_plasma, long_plasma$Accession %in% in_all_3)
gt3$Gene <- str_extract(gt3$Description, "(?<=GN=)[^ ]+")
gt3 <- gt3 %>% dplyr::select("Accession", "Gene", "Description")
gt3 <- unique(gt3)
gt3_table <- gt(gt3)
gt3_table
#save gt table with ggsave
#library(ragg)
#library(grid)
#ggsave("plots/3_in_all_table.svg", plot = gt3_table)

#gtsave(gt3_table, "plots/3_in_all_table.png")

# brain and plasma shared
bp14 <- subset(long_plasma, long_plasma$Accession %in% bp)
bp14$Gene <- str_extract(bp14$Description, "(?<=GN=)[^ ]+")
bp14 <- bp14 %>% dplyr::select("Accession", "Gene", "Description")
bp14 <- unique(bp14)
bp14_table <- gt(bp14)
bp14_table


#### Upset Plot ####
{
  brain <- read_excel("manuscript_data/NEW_AD_proteomics_R_results_Brain(Sept_2024).xlsx", sheet = "sig")
  csf <- read_excel("manuscript_data/NEW_AD_proteomics_R_results_CSF(Sept_2024).xlsx", sheet = "sig")
  plasma <- read_excel("manuscript_data/NEW_AD_proteomics_R_results_Plasma(Sept_2024).xlsx", sheet = "sig")
  brain_female <- read_excel("manuscript_data/AD_proteomics_R_results_Brain(Sept_2024).xlsx", sheet = "Female")
  brain_male <- read_excel("manuscript_data/AD_proteomics_R_results_Brain(Sept_2024).xlsx", sheet = "Male")
  CSF_female <- read_excel("manuscript_data/AD_proteomics_R_results_CSF(Sept_2024).xlsx", sheet = "Female")
  CSF_male <- read_excel("manuscript_data/AD_proteomics_R_results_CSF(Sept_2024).xlsx", sheet = "Male")
  plasma_female <- read_excel("manuscript_data/AD_proteomics_R_results_Plasma(Sept_2024).xlsx", sheet = "Female")
  plasma_male <- read_excel("manuscript_data/AD_proteomics_R_results_Plasma(Sept_2024).xlsx", sheet = "Male")
  
  #brain_female <- subset(brain_female, brain_female$F_AD_C_pval <0.05)
  #brain_male <- subset(brain_male, brain_male$M_AD_C_pval <0.05)
  #CSF_female <- subset(CSF_female, CSF_female$F_AD_C_pval <0.05)
  #CSF_male <- subset(CSF_male, CSF_male$M_AD_C_pval <0.05)
  #plasma_female <- subset(plasma_female, plasma_female$F_AD_C_pval <0.05)
  #plasma_male <- subset(plasma_male, plasma_male$M_AD_C_pval <0.05)
  
  proteins <- list(
    `Brain All` = c(brain$Accession),
    `Brain Female` = c(brain_female$Accession),
    `Brain Male` = c(brain_male$Accession),
    `CSF All` = c(csf$Accession), 
    `CSF Female` = c(CSF_female$Accession),
    `CSF Male` = c(CSF_male$Accession), 
    `Plasma All` = c(plasma$Accession),
    `Plasma Female` = c(plasma_female$Accession),
    `Plasma Male` = c(plasma_male$Accession)
  )
  
  prot_matrix <- list_to_matrix(proteins)
  colnames(prot_matrix) <- c("Brain", "CSF", "Plasma")
  prot_com_matrix <- make_comb_mat(prot_matrix, mode = "distinct")
  prot_com_matrix
  
  upset(prot_com_matrix)
  ?make_comb_mat
  m = make_comb_mat(proteins)
  UpSet(m, )
  ?UpSet
  
  #setting colors
  #this can also be done with hexadecimal
  main_bar_col <- c("violetred4")
  sets_bar_col <- c("turquoise4")
  matrix_col <- c("slateblue4")
  shade_col <- c("wheat4")
  
  
  as.m
  set_vars <- c("Brain", "CSF", "Plasma")
  pm2 <- as.data.frame(prot_matrix)
  pm3 <- t(pm2)
  upset(pm2, pm2[,1:3], 
        sets = set_vars,
        mb.ratio = mb_ratio1, 
        mainbar.y.label = "Counts by Pattern of Conditions", 
        sets.x.label = "Counts by Condition",
        order.by = "freq",
        show.numbers = T,
        point.size = 2, 
        line.size = 1,
        #text.scale=text_scale_options3,
        #main.bar.color = main_bar_col,
        #sets.bar.color = sets_bar_col,
        #matrix.color = matrix_col,
        #shade.color = shade_col)
  )
  
  x <- upset(fromList(proteins), 
             sets = colnames(proteins), nsets = 9, matrix.color = "gray25",
             main.bar.color = "steelblue", sets.bar.color = "wheat3", point.size = 2.5, 
             line.size = 1, text.scale = 2, set_size.show = T)
  ggsave("AD_upset.svg", height = 7, width = 7, units = "in")
}

#### Gender Venn ####
{
  
  sig <- read_excel("manuscript_data/AD_proteomics_R_results_Plasma(July_2024).xlsx", sheet = "sig")
  male <- read_excel("manuscript_data/AD_proteomics_R_results_Plasma(July_2024).xlsx", sheet = "Male")
  female <- read_excel("manuscript_data/AD_proteomics_R_results_Plasma(July_2024).xlsx", sheet = "Female")
  
  #Plasma
  a <- unlist(sig[,2])
  b <- unlist(male[,2])
  c <- unlist(female[,2])
  
  venn <- list(`AD vs C` = a, 
               `Male` = b, 
               `Female` = c)
  
  ggvenn(venn, show_percentage = F, text_size = 7, set_name_size = 10, 
         show_outside = "auto", fill_color = c("yellow3","steelblue","tomato3"), 
         fill_alpha = 0.6, stroke_size = 1.5)
  
  #gold2, tomato, skyblue
  #ggsave("plots/AD_sig_venn_proteins_gender.svg", width = 7.5, height = 7.5, units = "in")
  
  v.table <- venn(venn)
  sec <- attributes(v.table)
  
  #extract sections
  male_prots <- sec$intersections$Male
  female_prots <- sec$intersections$Female
}

#### DAVID ####
{
  ##### MF Brain/Plasma 15 #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\bp_15_molecular_function.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_MF_DIRECT = strsplit(as.character(GOTERM_MF_DIRECT), ",")) %>%
    unnest(GOTERM_MF_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_MF_DIRECT = gsub(".*~", "", GOTERM_MF_DIRECT)) %>%
    count(GOTERM_MF_DIRECT, name = "count") %>%
    filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_MF_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "goldenrod2") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Molecular Function GO Analysis of Statistically Signficant Proteins\nShared Between Brain and Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0))
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 5, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/mf_15_v3.svg", plot =  r, height = 4, width = 16.5)
  
  ##### BP Brain/Plasma 15  #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\bp_15_biological_process.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ",")) %>%
    unnest(GOTERM_BP_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_BP_DIRECT = gsub(".*~", "", GOTERM_BP_DIRECT)) %>%
    count(GOTERM_BP_DIRECT, name = "count") %>%
    filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_BP_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Biological Process GO Analysis of Statistically Signficant Proteins\nShared Between Brain and Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0))
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 1.1, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/bp_15_v3.svg", plot =  r, height = 6, width = 15)
  
  
  
  ##### MF 3 in all #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\mf_3inall.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_MF_DIRECT = strsplit(as.character(GOTERM_MF_DIRECT), ",")) %>%
    unnest(GOTERM_MF_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_MF_DIRECT = gsub(".*~", "", GOTERM_MF_DIRECT)) %>%
    count(GOTERM_MF_DIRECT, name = "count") #%>%
  #filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_MF_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "goldenrod2") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Molecular Function GO Analysis of Statistically Signficant Proteins\nIn Brain, CSF, and Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0))
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 1.1, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/mf_3inall_v3.svg", plot =  r, height = 10, width = 21)
  
  
  
  ##### BP 3 in all #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\bp_3inall.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ",")) %>%
    unnest(GOTERM_BP_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_BP_DIRECT = gsub(".*~", "", GOTERM_BP_DIRECT)) %>%
    count(GOTERM_BP_DIRECT, name = "count") #%>%
  #filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_BP_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Biological Process GO Analysis of Statistically Signficant Proteins\nIn Brain, CSF, and Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0))
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 1.1, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/bp_3inall_v3.svg", plot =  r, height = 15, width = 20)
  
  
  ##### MF Male #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\mf_male.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_MF_DIRECT = strsplit(as.character(GOTERM_MF_DIRECT), ",")) %>%
    unnest(GOTERM_MF_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_MF_DIRECT = gsub(".*~", "", GOTERM_MF_DIRECT)) %>%
    count(GOTERM_MF_DIRECT, name = "count") %>%
    filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_MF_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "goldenrod2") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Molecular Function GO Analysis of Statistically Signficant Proteins\nUnique to Male Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0))
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 1.1, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/mf_male_v3.svg", plot =  r, height = 6, width = 16)
  
  ##### BP Male #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\bp_male.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ",")) %>%
    unnest(GOTERM_BP_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_BP_DIRECT = gsub(".*~", "", GOTERM_BP_DIRECT)) %>%
    count(GOTERM_BP_DIRECT, name = "count") %>%
    filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_BP_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Biological Process GO Analysis of Statistically Signficant Proteins\nUnique to Male Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0))
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 1.1, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/bp_male_v3.svg", plot =  r, height = 8, width = 24)
  
  
  ##### MF Female #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\mf_female.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_MF_DIRECT = strsplit(as.character(GOTERM_MF_DIRECT), ",")) %>%
    unnest(GOTERM_MF_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_MF_DIRECT = gsub(".*~", "", GOTERM_MF_DIRECT)) %>%
    count(GOTERM_MF_DIRECT, name = "count") %>%
    filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_MF_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "goldenrod2") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Molecular Function GO Analysis of Statistically Signficant Proteins\nUnique to Female Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0)) #+
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 1.1, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/mf_female_v3.svg", plot =  r, height = 4, width = 16)
  
  ##### BP Female #####
  # Load the data
  data <- read.table(r"(C:\Users\andyb\Downloads\bp_female.txt)", header = TRUE, sep = "\t")
  
  # Separate the GOTERM_BP_DIRECT into unique IDs and count the hit
  df_split <- data %>%
    mutate(GOTERM_BP_DIRECT = strsplit(as.character(GOTERM_BP_DIRECT), ",")) %>%
    unnest(GOTERM_BP_DIRECT) %>%
    #remove text before ~ in the GO term
    mutate(GOTERM_BP_DIRECT = gsub(".*~", "", GOTERM_BP_DIRECT)) %>%
    count(GOTERM_BP_DIRECT, name = "count") %>%
    filter(count > 1)
  
  # Create the bar plot
  r <- ggplot(df_split, aes(x = GOTERM_BP_DIRECT, y = count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    theme_classic2() +
    coord_flip() +
    labs(title = "Biological Process GO Analysis of Statistically Signficant Proteins\nUnique to Female Plasma Samples", 
         x = "",
         y = "") +
    theme(axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 18, face = "bold"))+
    theme(plot.title = element_text(size = 24, face = "bold")) +
    scale_y_continuous(breaks = seq(floor(min(df_split$count-2)), ceiling(max(df_split$count)), by = 1)) +
    scale_y_continuous(expand = c(0, 0))
  #geom_text(aes(label = GOTERM_MF_DIRECT), hjust = 1.1, fontface = "bold", color = "dodgerblue4", size = 6)
  
  #r 
  r
  ggsave("manuscript_data/plots/DAVID_2/bp_female_v3.svg", plot =  r, height = 4, width = 15)
}












