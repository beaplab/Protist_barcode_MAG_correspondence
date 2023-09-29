## Script for processing raw output of proportionality calculation of shuffling. 

## Originally run under Linux 5.15.0-46-generic, Ubuntu 22.04.1 with x86-64 architecture, on a machine with 32 CPU and 125Gb RAM, with R version 4.1.2 (2021-11-01).  ("nisaba")

##*a comment below is on working paths in the local PCs*
## nisaba wd = , nisaba input wd = , nisaba output wd = , dazhbog wd = , dazhbog output wds = 


##inputs required:
# - final_taxa_ext_genus.csv
# - *rand_all_scores_relabund_inv_ext_shuff_scratch.tsv


##outputs produced:
# - summary_shuffle_stats_all.tsv - summary stats for all iterations of shufflings; input for rho_calc_plots.Rmd
# - se_summary_shuffle_stats_all.tsv - statistics for SD accumulation plot; input for se_shuffle_plots.Rmd



##Library the packages (some of them not actually required here but needed for producing plots; those are commented)

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


#principle of naming files was as follows:
#paste0("/scratch/data1/dzavadska/shuffle/",metric, i,  sg, iteration, "_rand_","all_scores_relabund_inv_ext_shuff_scratch.tsv")
#
setwd("/scratch/data1/dzavadska/shuffle/")


##############################################################################################
#======================================BLOCK 1===============================================#
##############################################################################################
#PRODUCING THE PRIMARY BIG DF WITH ALL THE DATA


best3 <- list.files(pattern = "rand_all_scores_relabund_inv_ext_shuff_scratch.tsv")

#empty files (of size 16K bcuz of column names) generated in the past in directory, they are normally deleted, however, to avoid inconveniences, let's filter by file size.
best3 <- best3[sapply(best3, file.size) > 16]

#subsetting to only 2 different metrics that are included in main graph not to overload everything; UPD====> in this version of script run on nisaba this step is not required
#best3 <- c(best3[grep("rho", best3)], best3[grep("cor", sub("scor","",best3))])
best3 <- c(best3[grep("rho", best3)])
#this is to save computational time!!! UPD====> in this version of script run on nisaba this step is not required
#best3 <- c(best3[grep("5", best3)])

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
#backup0
#df_backup <- df #UPD====> in this version of script run on nisaba this step is not required

df$ftop <- as.factor(df$tops) 
df <- df[which(df$tops<=10),]

df$organism <- str_replace(df$organism, "_rand_all_scores_relabund_inv_ext_shuff_scratch.tsv", "")
#backup1
df$taxlevel <- df$organism

df$organism <- str_replace(df$organism, "supergroup|taxogroup_1|taxogroup_2|taxogroup_middle|genus", "")
#backup2

#workaround for New_***02|01 groups
#levels(as.factor(df$organism[grep("_01|_02", df$organism)]))
#levels(as.factor(df$organism[grep("_02", df$organism)]))
groups2 <- grep("_02", df$organism)
groups1 <- grep("_01", df$organism)

df$iteration <- extract_numeric(df$organism)

#workaround for New_***02|01 groups
df$iteration[groups2] <-  str_replace(df$iteration[groups2], "2", "")
df$iteration[groups1] <-  str_replace(df$iteration[groups1], "1", "")

df$organism <- str_replace(df$organism, as.character(df$iteration), "")
#
#grep("supergroup|taxogroup_1|taxogroup_2|taxogroup_middle|genus", df$organism) #checks if taxolevels were successfully removed
#2 lines below are workaround for Pseudo-nitschia bitch hyphen
df$iteration <- str_replace(df$iteration, "-", "")
df$organism <- str_replace(df$organism, as.character(df$iteration), "")
#workaround for New_MAST-4 bitch hyphen
df$organism[which(df$organism=="rhoNew_MAST")] <- "rhoNew_MAST-4"
df$organism[which(df$organism=="corNew_MAST")] <- "corNew_MAST-4"
df$organism[which(df$organism=="phiNew_MAST")] <- "phiNew_MAST-4"
df$organism[which(df$organism=="phsNew_MAST")] <- "phsNew_MAST-4"
df$organism[which(df$organism=="vlrNew_MAST")] <- "vlrNew_MAST-4"



df$taxlevel <- str_replace(df$taxlevel, df$organism, "")
df$taxlevel <- str_replace(df$taxlevel, as.character(df$iteration), "")

#check if number replacements didnt affect taxogroup_1/2 numbering --> its ok
#df[intersect(which(df$iteration==1),which(df$taxlevel=="taxogroup_1")),]



vocabulary <- fread("final_taxa_ext_genus.csv")


df$metrics<-NA
df$metrics <- substr(df$organism, 1, 3)  

diffmetr <- c("phi","rho","phs","cor","vlr")
df <- df[df$metrics%in%diffmetr,]
#levels(as.factor(df$metrics)) check if all metrics was included

df$organism <- sub("^...","",df$organism)

######################OPTIONALLY - if you want
#df$taxogroup<-NA

#for (i in c(1:nrow(df))){
#  colnumer<-which(colnames(vocabulary)==df$taxlevel[i])
#  vocabulary_dummy <- vocabulary %>% select(all_of(colnumer))
#  vocabulary_dummy$taxa_SMAGs <- vocabulary$taxa_SMAGs
#  df$taxogroup[i]<-as.character(vocabulary_dummy[which(vocabulary$taxa_SMAGs==df$organism[i]),1][1])
#}
#backup3

#adding some non-essential variables to the table UPD====>  columns with this variable do not seem to be required by now in this dataframe. NOTE====> Keep it just in case by now, DELETE IN THE FUTURE MAYBE
#df$fSMAG <- paste0(df$SMAG, df$metrics)
#df$taxlfSMAG <- paste0(df$SMAG, df$metrics, df$taxlevel)

#workaround for New_MAST-4 bitch hyphen continues
df$iteration[which(df$organism=="New_MAST-4")] <- str_replace(df$iteration[which(df$organism=="New_MAST-4")], "4", "")
df$taxlevel[which(df$organism=="New_MAST-4")] <-str_replace(df$taxlevel[which(df$organism=="New_MAST-4")], as.character(df$iteration[which(df$organism=="New_MAST-4")]), "")




##############################################################################################
#======================================BLOCK 2===============================================#
##############################################################################################
#THIS BLOCK TAKES ONLY RECORDS FOR 2 LOWEST TAXONOMIC LEVELS (BY V9 TAXONOMY)
# NOTE====> THIS BLOCK MIGHT NEED A BIT OF MANUAL MODIFICATIONIF FOR SOME GROUPS SOME SPECIFIC TAXLEVELS SHOULD BE PLOTTED



###############with lowest taxonomic rank up to taxogroup_2
df_no_super <- df[which(df$taxlevel != "supergroup"),]
df_no_tax1 <- df_no_super[which(df_no_super$taxlevel != "taxogroup_1"),]


df_na_vocab <- data.frame()


for (i in c(4:(ncol(vocabulary)))){
  #print(i)
  df_low_taxlevel_rank <- c()
  df_low_taxlevel_group <- c()
  
  
  vocabulary_dummy <- vocabulary %>% select(all_of(i))
  df_low_taxlevel_group <- c(df_low_taxlevel_group,vocabulary$taxa_SMAGs[which(is.na(vocabulary_dummy[,1]==TRUE))])
  df_low_taxlevel_rank <- c(df_low_taxlevel_rank, colnames(vocabulary)[i-1])
  
  taxlevel <- rep(df_low_taxlevel_rank,length(vocabulary$taxa_SMAGs[which(is.na(vocabulary_dummy[,1]==TRUE))]))
  data <- data.frame(matrix(nrow=length(taxlevel)))
  
  data$group <- df_low_taxlevel_group
  data$taxlevel <- taxlevel
  
  df_na_vocab<- rbind.fill(df_na_vocab, data)
  
}
#df_na_vocab


group <- c()
taxlevel <- c()
for(i in unique(df_na_vocab$group)){
  #print(i)
  #print(df_na_vocab$taxlevel[which(df_na_vocab$group == i)][1])
  taxlevel <- c(taxlevel, df_na_vocab$taxlevel[which(df_na_vocab$group == i)][1])
  group <- c(group, i)
}

df_na_fin <- data.frame(matrix(nrow=length(taxlevel)))
df_na_fin$group <- group
df_na_fin$taxlevel <- taxlevel



subs_trunc_plus <- data.frame()

ranks <- c("supergroup", "taxogroup_1", "taxogroup_2", "taxogroup_middle", "genus")

for(i in c(1:length(df_na_fin$group))){
  #print(i)
  #print(df_na_fin$group[i])
  subs <- df[which(df$organism == df_na_fin$group[i]),]
  subs_lvl1 <- subs[which(subs$taxlevel == df_na_fin$taxlevel[i]),]
  
  if (df_na_fin$taxlevel[i]!="supergroup") {
    higher_rank <- print(ranks[which(ranks == df_na_fin$taxlevel[i])-1])
    subs_lvl2 <- subs[which(subs$taxlevel == higher_rank),]
    subs_trunc_plus <- rbind.fill(subs_trunc_plus,subs_lvl1,subs_lvl2)
  } else
  {subs_trunc_plus <- rbind.fill(subs_trunc_plus,subs_lvl1)}
  
  
}
#check if it worked; should be output "taxogroup_2" AND "taxogroup_1"
#unique(subs_trunc_plus$taxlevel[which(subs_trunc_plus$organism=="New_MAST-4")])


#as side product, df_na_fin was written out. This will be used for SMAGs scatterplot input
#fwrite(df_na_fin[,2:3], "df_na_fin_for_scatterplots.tsv")

df_lvl1 <- df[-which(df$organism %in% df_na_fin$group),]
df_lvl2 <- df_lvl1[which(df_lvl1$taxlevel == "genus"),]
df_lvl3 <- df_lvl1[which(df_lvl1$taxlevel == "taxogroup_middle"),]

to_plot <- rbind.fill(subs_trunc_plus,df_lvl2,df_lvl3)


#to_plot_rho <- to_plot[which(to_plot$metrics == "rho"),]
#to_plot_cor <- to_plot[which(to_plot$metrics == "cor"),]
#to_plot_rho <- to_plot[which(to_plot$metrics == "phi"),]
#to_plot_cor <- to_plot[which(to_plot$metrics == "phs"),]
#to_plot_cor <- to_plot[which(to_plot$metrics == "vlr"),]

#####################################################################################################
###############OPTION FOR QUICK PLOTS- instead of assigning taxogroup to everything, we can assign it only to those taxlevels that are plotted. More generally, we can even skip the final taxogroup assignation, as in the end, final plots take the taxlevel as am argument, not taxogroup. So uncomment if you need it. 
#####################################################################################################

#df$taxogroup<-NA

#for (i in c(1:nrow(df))){
#  colnumer<-which(colnames(vocabulary)==df$taxlevel[i])
#  vocabulary_dummy <- vocabulary %>% select(all_of(colnumer))
#  vocabulary_dummy$taxa_SMAGs <- vocabulary$taxa_SMAGs
#  df$taxogroup[i]<-as.character(vocabulary_dummy[which(vocabulary$taxa_SMAGs==df$organism[i]),1][1])
#}

#####################################################################################################
###############OPTION FOR QUICK PLOTS- instead of assigning taxogroup to everything, we can assign it only to those taxlevels that are plotted. More generally, we can even skip the final taxogroup assignation, as in the end, final plots take the taxlevel as am argument, not taxogroup. So uncomment if you need it. 
#####################################################################################################


ranks <- c("supergroup", "taxogroup_1", "taxogroup_2", "taxogroup_middle", "genus")

################for big to_plot dataframe with all the metrics

#to_plot$ittaxlfSMAG <- paste0(to_plot$iteration, to_plot$taxlfSMAG)

organism <- c()
taxlevel <- c()
Min <- c()
Max <- c()
Mean <- c()
Se <- c()
top <- c()
metrics <- c()
iteration_count <- c()


to_plot$grouping <- paste0(to_plot$organism, to_plot$taxlevel, to_plot$metrics, to_plot$tops)

##removing non-distinct rows from the to_plot; because it results in many non-unique rows for iterations and md5sum doesnt matter in the end because it is shuffling
to_plot <- distinct(to_plot)

for (X in levels(as.factor(to_plot$grouping))) {
  rho <- to_plot$rho[which(to_plot$grouping == X)]
  Min <- c(Min, min(rho))
  Max <- c(Max, max(rho))
  Mean <- c(Mean, mean(rho))
  Se <- c(Se, sd(rho))
  
  organism <- c(organism, unique(to_plot$organism[which(to_plot$grouping == X)]))
  taxlevel <- c(taxlevel, unique(to_plot$taxlevel[which(to_plot$grouping == X)]))
  top <- c(top, unique(to_plot$top[which(to_plot$grouping == X)]))
  metrics <- c(metrics, unique(to_plot$metrics[which(to_plot$grouping == X)]))
  
  iteration_count <- c(iteration_count, length(unique(to_plot$iteration[which(to_plot$grouping == X)])))
}


summary_shuffle_stats_all <- data.frame(
  "organism" = organism, 
  "taxlevel"= taxlevel,
  "Min"=Min,
  "Max"=Max,
  "Mean"=Mean,
  "Se"=Se,
  "top" = top,
  "metrics" = metrics,
  "iteration_count" = iteration_count
)

#usually duplicate rows do not appear at this point, but just in case
summary_shuffle_stats_all <- distinct(summary_shuffle_stats_all)


fwrite(summary_shuffle_stats_all, "summary_shuffle_stats_all.tsv")




##############################################################################################
#======================================BLOCK 3===============================================#
##############################################################################################
#THIS BLOCK PREPARES DATA FOR SE PLOT

organism_se <- c()
taxlevel_se <- c()
Min_se <- c()
Max_se <- c()
Mean_se <- c()
Se_se <- c()
top_se <- c()
metrics_se <- c()
iteration_count_se <- c()

itr <-c()

to_plot_backup <- to_plot

to_plot <- to_plot[which(to_plot$top == 1),] 


  
for (X in levels(as.factor(to_plot$grouping))) {
  for (XX in c(1:length(levels(as.factor(to_plot$iteration))))) {
  itr <- c(1:XX)
  
  rho <- to_plot$rho[intersect(which(to_plot$iteration %in% itr),which(to_plot$grouping == X))]
  if (length(rho)==0) {print(XX)}
  
  Min_se <- c(Min_se, min(rho))
  Max_se <- c(Max_se, max(rho))
  Mean_se <- c(Mean_se, mean(rho))
  Se_se <- c(Se_se, sd(rho))
  
  organism_se <- c(organism_se, unique(to_plot$organism[which(to_plot$grouping == X)]))
  taxlevel_se <- c(taxlevel_se, unique(to_plot$taxlevel[which(to_plot$grouping == X)]))
  top_se <- c(top_se, unique(to_plot$top[which(to_plot$grouping == X)]))
  metrics_se <- c(metrics_se, unique(to_plot$metrics[which(to_plot$grouping == X)]))
  
  iteration_count_se <- c(iteration_count_se, XX)
  
  }
}

summary_shuffle_stats_all_se <- data.frame(
  "organism" = organism_se, 
  "taxlevel"= taxlevel_se,
  "Min"=Min_se,
  "Max"=Max_se,
  "Mean"=Mean_se,
  "Se"=Se_se,
  "top" = top_se,
  "metrics" = metrics_se,
  "iteration_count" = iteration_count_se
)

#usually duplicate rows do not appear at this point, but just in case
summary_shuffle_stats_all_se <- distinct(summary_shuffle_stats_all_se)


fwrite(summary_shuffle_stats_all_se, "se_summary_shuffle_stats_all.tsv")



