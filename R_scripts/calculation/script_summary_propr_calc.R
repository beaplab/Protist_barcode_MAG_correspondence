
## Script for processing raw output of proportionality calculation of shuffling. 

## Originally run under Linux 5.15.0-46-generic, Ubuntu 22.04.1 with x86-64 architecture, on a machine with 32 CPU and 125Gb RAM, with R version 4.1.2 (2021-11-01).  ("nisaba")

##*a comment below is on working paths in the local PCs*
## nisaba wd = , nisaba input wd = , nisaba output wd = , dazhbog wd = , dazhbog output wds = 


##inputs required:
# - final_taxa_ext_genus.csv
# - *all_scores_relabund_upd_ext_gen.tsv

#outputs produced:
# - summary_propr_calc_stats_all.tsv


library(data.table)
#library(ggplot2)
library(reshape2)
library(tidyverse)
library(plyr)
library(dplyr)
#library(microshades)
#library(scales)
#library(ggbreak) 
#library(patchwork)
library(ggpubr)
#library(gplots)

setwd("/home/dzavadska/preprocess_newv9/")

best3 <- list.files(pattern = "all_scores_relabund_upd_ext_gen.tsv")

#subsetting to only 2 different metrics that are included in main graph not to overload everything
#best3 <- c(best3[grep("rho", best3)], best3[grep("cor", sub("scor","",best3))])

df <-data.frame()



for (i in c(1:length(best3))) {
  data <- fread(best3[i])
  
  data <- split(data,data[,SMAG]) %>%
    lapply(function(X){
      X[order(rho,decreasing = T)]
    }) %>% rbindlist() 
  
  x <- seq(1,length(unique(data$md5sum)),1)
  tops <- rep(x, nrow(data)/length(unique(data$md5sum)))
  
  data$tops <- tops
  data$organism <- rep(best3[i], nrow(data))
  
  df<- rbind.fill(df, data)
  
}

df$organism <- str_replace(df$organism, "all_scores_relabund_upd_ext_gen.tsv", "")

df$taxlevel <- df$organism

df$organism <- str_replace(df$organism, "supergroup|taxogroup_1|taxogroup_2|taxogroup_middle|genus", "")
#grep("supergroup|taxogroup_1|taxogroup_2|taxogroup_middle|genus", df$organism) #checks if taxolevels were successfully removed

df$taxlevel <- str_replace(df$taxlevel, df$organism, "")

df$ftop <- as.factor(df$tops) 
df <- df[which(df$tops<=10),]

vocabulary <- fread("final_taxa_ext_genus.csv")


df$metrics<-NA

df$metrics <- substr(df$organism, 1, 3)  


diffmetr <- c("phi","rho","phs","cor","vlr")
df <- df[df$metrics%in%diffmetr,]

levels(as.factor(df$metrics))

df$organism <- sub("^...","",df$organism)


df$taxogroup<-NA

#for (s in colnames(vocabulary)[3:7]){
# colnumer<-which(colnames(vocabulary)[3:7]==df$taxlevel)

for (i in c(1:nrow(df))){
  colnumer<-which(colnames(vocabulary)==df$taxlevel[i])
  vocabulary_dummy <- vocabulary %>% select(all_of(colnumer))
  vocabulary_dummy$taxa_SMAGs <- vocabulary$taxa_SMAGs
  df$taxogroup[i]<-as.character(vocabulary_dummy[which(vocabulary$taxa_SMAGs==df$organism[i]),1][1])
}
#}

df$fSMAG <- paste0(df$SMAG, df$metrics)
df$taxlfSMAG <- paste0(df$SMAG, df$metrics, df$taxlevel)

fwrite(df, "summary_propr_calc_stats_all.tsv")