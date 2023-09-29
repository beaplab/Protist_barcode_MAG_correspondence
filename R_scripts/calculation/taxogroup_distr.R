## Script for creating a summary stats on number of SMAGs and V9s included in a dataset. 
#Later those numbers are input for "final_simulation_clean_pipeline" script

## Originally run under Linux 5.15.0-46-generic, Ubuntu 22.04.1 with x86-64 architecture, on a machine with 32 CPU and 125Gb RAM, with R version 4.1.2 (2021-11-01).  ("nisaba")

##*a comment below is on working paths in the local PCs*
## nisaba wd = setwd("/home/dzavadska/preprocess_newv9"), nisaba input wd = setwd("/home/dzavadska/preprocess_newv9"), nisaba output wd = setwd("/home/dzavadska/preprocess_newv9"), dazhbog wd = , dazhbog output wds = 


##inputs required:
# - 18S_V9_1476_samples.OTU.filtered.cleaved.nosubstringOTUs.taxo
# - Table_S04_distributions_nr_SMAGs_for_SAGs_for_all.xlsx !!!IMPORTANT  - THE UPDATED TABLE WITH "MOCK" SMAGS CONTROL GENRES

#outputs produced:
# - taxa_counts.tsv --- used as input for simulation pipeline propr calculation; stored in  data_input AND data_output


##upload packages
library(stringr)
library(openxlsx)
library(data.table)
library(magrittr)
library(GGally)
library(propr)
##upload packages
################################
## set working directory
#for nisaba
setwd("/home/dzavadska/preprocess_newv9")

#for dazhbog
#

## set working directory


####The piece of data import; is essentially similar to the one in proportionality calculation script
#V9 new taxonomy assignment data
datataxa<- fread("18S_V9_1476_samples.OTU.filtered.cleaved.nosubstringOTUs.taxo")

taxa_SMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs_for_SAGs_for_all.xlsx", sheet = 4, startRow = 1, colNames = TRUE) %>% data.frame
#taxa_SMAGs - contains info on which SMAG belongs to which taxogroup

tmp <- data.frame(levels(factor(taxa_SMAGs[,2]))) 
#subset taxogroup colunm and assign factor to it

taxa_SMAGs_melted <- merge(taxa_SMAGs,tmp,by.x="AA_Best_taxonomy_GENRE",by.y="levels.factor.taxa_SMAGs...2...",all.x = T)
taxa_SMAGs_melted <- taxa_SMAGs_melted[,c("AA_Best_taxonomy_GENRE", "AA_SMAG")]
#Now we have a version of taxa SMAGs where each row contains a SMAG ID and a taxogroup where it is assigned

vocabulary <- fread("final_taxa_ext_genus.csv",sep=",")
taxa_list <- c(levels(factor(taxa_SMAGs_melted[,1]))[2:length(levels(factor(taxa_SMAGs_melted[,1])))])
list_of_all_old <- (c(taxa_list, c(unique(vocabulary$taxa_SMAGs))))
list_of_all_old <- list_of_all_old[duplicated(c(taxa_list, c(unique(vocabulary$taxa_SMAGs))))==T]


########### calculation for SMAGs

dataSMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs_for_SAGs_for_all.xlsx", sheet = 7, startRow = 1, colNames = TRUE) %>% data.table
dataSMAGs <- dataSMAGs[1:713,]
dataV9_melted1 <- fread("dataV9_meltedLOG-TRANSFORMED.tsv")


SMAGs_taxo_distr <- c()
SMAGs_genre <- c()
#an alternative would be to create another loop by uncommenting the line below to automatically iterate over taxlevels, but found it faster to do it manually. (moreover, it does not matter for SMAGs)
#for (d in c(3:length(colnames(vocabulary)))){
d <- 3
list_of_all <- c()
sg <- colnames(vocabulary)[d]
for (o in list_of_all_old ) {if(is.na(vocabulary[which(vocabulary$taxa_SMAGs == o), which(colnames(vocabulary)==colnames(vocabulary)[d])])!=TRUE){list_of_all <- c(list_of_all, o)}}
#i <- "Attheya"

for(i in list_of_all[1:(length(list_of_all))]) {
  #subsetting SMAGs dataset by the taxogroup
  
  tokeep <- c(taxa_SMAGs_melted[which(taxa_SMAGs_melted$AA_Best_taxonomy_GENRE == i),2])
  
  SMAGs_taxo_distr <- c(SMAGs_taxo_distr, length(tokeep))
  SMAGs_genre <- c(SMAGs_genre, i)
}    
#}

########### calculation for SMAGs

#a dataframe that will store stats
SMAGs_number_of_g <- data.frame("GENRE" = SMAGs_genre, "Number_of_SMAGs" = SMAGs_taxo_distr, "supergroup" = NA, "taxogroup_1" = NA, "taxogroup_2" = NA, "taxogroup_middle" = NA, "genus" = NA)


########### calculation for V9s
dataV9 <- fread("18S_V9_1476_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.dereplicated.table (3)")
dataV9 <- dataV9[,c(1, 5:1088)] 
  

for (d in c(3:length(colnames(vocabulary)))){
#uncomment line below to run test xmple
#d <- 7
list_of_all <- c()
sg <- colnames(vocabulary)[d]
V9_taxo_distr <- c()

for (o in list_of_all_old ) {if(is.na(vocabulary[which(vocabulary$taxa_SMAGs == o), which(colnames(vocabulary)==colnames(vocabulary)[d])])!=TRUE){list_of_all <- c(list_of_all, o)}}

for(i in list_of_all[1:(length(list_of_all))]) {
  #uncomment line below to run test xmple
  #i <- "Attheya"
  taxogroup_i <- data.frame(vocabulary)[which(vocabulary == i)[1], d]
  amplicons <- datataxa[grep(taxogroup_i, datataxa$taxonomy),"amplicon"]

  if (length(which(dataV9$amplicon %in% amplicons$amplicon))==0) {
    V9_taxo_distr <- c(V9_taxo_distr, 0)
  } 
  else {
  dataV9_tmp <- dataV9[which(dataV9$amplicon %in% amplicons$amplicon),]
  V9_taxo_distr <- c(V9_taxo_distr, nrow(dataV9_tmp))
  
  }
  
}    

SMAGs_number_of_g[,d] <- V9_taxo_distr

}

########### calculation for V9s

fwrite(SMAGs_number_of_g, "taxa_counts.tsv", sep = "\t")



