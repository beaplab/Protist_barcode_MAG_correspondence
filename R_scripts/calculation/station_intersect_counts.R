

#barcode-SMAG co-occurences 

## WARNING. the original script was run on nder Linux 5.15.0-46-generic, Ubuntu 22.04.1 with x86-64 architecture, on a machine with 32 CPU and 125Gb RAM, with R version 4.1.2 (2021-11-01). This is a computationally demanding which is likely not supposed to run on the local computers.

#datataxa[grep("Attheya", datataxa$taxonomy),"amplicon"]

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

setwd("/home/dzavadska/preprocess_newv9")

## set working directory

datataxa<- fread("18S_V9_1476_samples.OTU.filtered.cleaved.nosubstringOTUs.taxo")

taxa_SMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs.xlsx", sheet = 4, startRow = 1, colNames = TRUE) %>% data.frame
#taxa_SMAGs - contains info on which SMAG belongs to which taxogroup

tmp <- data.frame(levels(factor(taxa_SMAGs[,2]))) 
#subset taxogroup colunm and assign factor to it

taxa_SMAGs_melted <- merge(taxa_SMAGs,tmp,by.x="AA_Best_taxonomy_GENRE",by.y="levels.factor.taxa_SMAGs...2...",all.x = T)
taxa_SMAGs_melted <- taxa_SMAGs_melted[,c("AA_Best_taxonomy_GENRE", "AA_SMAG")]
#Now we have a version of taxa SMAGs where each row contains a SMAG ID and a taxogroup where it is assigned



########### Creating log-transformed table of SMAGs counts

dataSMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs.xlsx", sheet = 7, startRow = 1, colNames = TRUE) %>% data.table
dataSMAGs <- dataSMAGs[1:713,]


dataSMAGs_melted <- melt(dataSMAGs,id.vars="SMAG")

context_SMAGs <- fread("TAGs-18S-V4_NAME-PROJ-MAT-READ-CAB_nico.list",header = F)
tmp <- fread("TARA_CONTEXT_95_SEQ_STATIONS.txt.gz")
tmp <- tmp[,list(Sample.ID,Sample.mat)]
tmp[,Sample.ID:=sub("TARA_","",Sample.ID)]
tmp <- merge(context_SMAGs,tmp,by.x="V5",by.y="Sample.ID",all.x=T)
tmp <- tmp[,list(V1,Sample.mat)]
tmp <- unique(tmp)
tmp <- tmp[!duplicated(V1)]
dataSMAGs_melted <- merge(dataSMAGs_melted,tmp,by.x="variable",by.y="V1",all.x = T)
dataSMAGs_df <- dcast(Sample.mat~SMAG,value.var="value",data=dataSMAGs_melted[!is.na(Sample.mat)],fun.aggregate = mean) %>%
  data.frame(row.names = "Sample.mat")






dataV9 <- fread("18S_V9_1476_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.dereplicated.table (3)")
dataV9 <- dataV9[,c(1, 5:1088)] 

#dataV9_melted[,value2:=value/sum(value),by=variable]
####################################################################

##upd:error fixed #not all amplicons present in datataxa are also present in dataV9. THis is exeptionally strange, and should be investigated in more detail later. I will treat it as artifact of analysis by now. but weird. the line below is a solution to treat such cases

#dataV9 <- dataV9[taxogroup==c(taxogroup_i)[1]&pid>90]

x <- apply(dataV9[,.SD,.SDcols=grep("TARA",colnames(dataV9),value=T)],1,function(X) sum(X>0))
dataV9 <- dataV9[x>100]


dataV9_melted <- dataV9[,.SD,.SDcols=c("amplicon",grep("TARA",colnames(dataV9),value=T))] %>% melt(id.vars="amplicon")

contextV9 <- fread("TARA_CONTEXT_95_SEQ_STATIONS_V9.tsv")
dataV9_melted <- merge(dataV9_melted,contextV9[,list(Pangaea_sample_id,Sample.mat)],by.x="variable",by.y="Pangaea_sample_id")
dataV9_melted[,variable:=NULL]
dataV9_df_relabund <- dcast(Sample.mat~amplicon,value.var="value",data=dataV9_melted,fun.aggregate = mean) %>%
  data.frame(row.names = "Sample.mat",check.names = F)


#dataV9_df_relabund <- dataV9_df_relabund[row.names(dataSMAGs_df),]

tmp1 <- dataSMAGs_df

##############PROPR PACKAGE EDITION STARTS HERE################
tmp1 <- na.omit(tmp1)
#deleting NAs; 0s will be automatically replaced by propr package function
#tmp1[is.na(tmp1)==TRUE] <- 0
tmp1 <- merge(as.data.frame(tmp1),as.data.frame(dataV9_df_relabund),by="row.names") %>%
  data.frame(row.names = "Row.names",check.names = F)



##preparing the loop to calculate station number, subsetting within the taxogroup -- to save calculation time.

SMAG_names <- colnames(tmp1)[grep("TARA|Metagenome", colnames(tmp1))]

V9_names <- colnames(tmp1)[-grep("TARA|Metagenome", colnames(tmp1))]

SMAGs_r <-c()
V9s_r <- c()
count_r <- c()



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


#count the number of intersections per supergroup. with few exceptions, supergroups contain all V9 OTUs that are also within the lower taxonomic ranks.
d <- 3

list_of_all <- c()
sg <- colnames(vocabulary)[d]

for (o in list_of_all_old ) {if(is.na(vocabulary[which(vocabulary$taxa_SMAGs == o), which(colnames(vocabulary)==colnames(vocabulary)[d])])!=TRUE){list_of_all <- c(list_of_all, o)}}

for(i in list_of_all[1:(length(list_of_all))]) {
  
  SMAG_names_subs <- c(taxa_SMAGs_melted[which(taxa_SMAGs_melted$AA_Best_taxonomy_GENRE == i),2])
taxogroup_i <- data.frame(vocabulary)[which(vocabulary == i)[1], d]
  
  
V9_names_subs <-  datataxa[grep(taxogroup_i, datataxa$taxonomy),"amplicon"]
  
  
for (y in SMAG_names_subs) {
  for (e in V9_names_subs$amplicon) {
    if (e %in% colnames(tmp1) == FALSE) {next   }
    
    SMAGs_r <-c(SMAGs_r, y)
    V9s_r <- c(V9s_r, e)
    count_r <- c(count_r, length(intersect(which(tmp1[,y]!=0), which(tmp1[,e]!=0) )))
    
  }
}

  
}

df_fin <- data.frame("SMAG" = SMAGs_r, "md5sum" = V9s_r, "stcount" = count_r)

fwrite(df_fin,"station_intersect_counts.tsv")
