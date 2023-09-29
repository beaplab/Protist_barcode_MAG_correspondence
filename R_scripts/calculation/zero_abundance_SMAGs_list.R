
library(dplyr)
library(openxlsx)
library(data.table) 

taxa_SMAGs <- read.xlsx("/home/dzavadska/Data/MAGs/MAGs_pre-final/Ext_tax_bug_corr_SMAGs_genus/Table_S04_distributions_nr_SMAGs_for_SAGs_for_all.xlsx", sheet = 7, startRow = 1, colNames = TRUE) %>% data.frame

taxa_SMAGs_taxonomy <- read.xlsx("/home/dzavadska/Data/MAGs/MAGs_pre-final/Ext_tax_bug_corr_SMAGs_genus/Table_S04_distributions_nr_SMAGs_for_SAGs_for_all.xlsx", sheet = 4, startRow = 1, colNames = TRUE) %>% data.frame


x <- apply(taxa_SMAGs[,.SD,.SDcols=grep("X",colnames(taxa_SMAGs),value=T)],1,function(X) sum(X>0))

x <- apply(taxa_SMAGs[,c(2:ncol(taxa_SMAGs))],1,function(X) sum(X>0))

length(which(x == 0))

vocab <- fread("/home/dzavadska/Data/MAGs/MAGs_pre-final/Ext_tax_bug_corr_SMAGs_genus/final_taxa_ext_genus.csv")

y <- c(vocab[,1])

zero <- which(x == 0)

taxa_SMAGs$SMAG[zero]



######################################################old version for another spreadsheet
protist_SMAGs <- taxa_SMAGs_taxonomy$AA_SMAG[which(as.vector(taxa_SMAGs_taxonomy[,2]) %in% y$taxa_SMAGs == TRUE)]

z <- zero[which(c(as.vector(taxa_SMAGs[which(x == 0),2]) %in% y$taxa_SMAGs) == TRUE)]

w <-  zero[which(is.na(as.vector(taxa_SMAGs[which(x == 0),2])) == TRUE)]
taxa_SMAGs[c(z,w),1]
 

c("TARA_AOS_82_MAG_00178", "TARA_ION_45_MAG_00153", "TARA_IOS_50_MAG_00155", "TARA_MED_95_MAG_00394", "TARA_MED_95_MAG_00399", "TARA_MED_95_MAG_00442", "TARA_MED_95_MAG_00470", "TARA_MED_95_MAG_00516", "Metagenome_centric_SAG_TOSAG00_7",  
   "Metagenome_centric_SAG_TOSAG00_8", "Metagenome_centric_SAG_TOSAG23_6", "Metagenome_centric_SAG_TOSAG23_9", "Metagenome_centric_SAG_TOSAG39_1", "Metagenome_centric_SAG_TOSAG39_4", "Metagenome_centric_SAG_TOSAG39_6", "Metagenome_centric_SAG_TOSAG41_1",
   "Metagenome_centric_SAG_TOSAG41_10", "Metagenome_centric_SAG_TOSAG41_11", "Metagenome_centric_SAG_TOSAG41_9", "Metagenome_centric_SAG_TOSAG46_2",  "Metagenome_centric_SAG_TOSAG47_5","Metagenome_centric_SAG_TOSAG47_7", "Metagenome_centric_SAG_TOSAG48_3",
   "TARA_AOS_82_MAG_00159","TARA_AOS_82_MAG_00170", "TARA_ARC_108_MAG_00231", "TARA_ARC_108_MAG_00242", "TARA_ARC_108_MAG_00271", "TARA_ARC_108_MAG_00300", "TARA_ION_45_MAG_00172", "TARA_ION_45_MAG_00176" ,"TARA_ION_45_MAG_00190", "TARA_ION_45_MAG_00192", 
   "TARA_ION_45_MAG_00195", "TARA_ION_45_MAG_00198", "TARA_IOS_50_MAG_00137", "TARA_IOS_50_MAG_00167", "TARA_MED_95_MAG_00393", "TARA_MED_95_MAG_00397", "TARA_MED_95_MAG_00416", "TARA_MED_95_MAG_00450", "TARA_MED_95_MAG_00467", "TARA_MED_95_MAG_00473", 
   "TARA_MED_95_MAG_00500", "TARA_MED_95_MAG_00505", "TARA_PSE_93_MAG_00261", "Metagenome_centric_SAG_TOSAG41_8" )