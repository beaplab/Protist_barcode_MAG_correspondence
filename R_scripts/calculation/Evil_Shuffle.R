## WARNING. the original script was run on nder Linux 5.15.0-46-generic, Ubuntu 22.04.1 with x86-64 architecture, on a machine with 32 CPU and 125Gb RAM, with R version 4.1.2 (2021-11-01). This is a computationally demanding which is likely not supposed to run on the local computers.

#datataxa[grep("Attheya", datataxa$taxonomy),"amplicon"]

##upload packages
library(stringr)
library(openxlsx)
library(data.table)
library(magrittr)
library(GGally)
library(propr)
library(picante)
library(Rfast)
##upload packages
################################
## set working directory

setwd("/home/dzavadska/preprocess_newv9/")


###iqlr transformation function from propr package
##
ivar2index <- function(counts){
  ct <- counts
  if(any(counts == 0)){
    message("Alert: Replacing 0s with next smallest value.")
    zeros <- counts == 0
    counts[zeros] <- min(counts[!zeros])
  }
  
  counts.clr <- apply(log(counts), 1, function(x){ x - mean(x) })##
  counts.var <- apply(counts.clr, 1, var)##
  quart <- stats::quantile(counts.var) # use features with unextreme variance
  use <- which(counts.var < quart[4] & counts.var > quart[2])
  
  
  #return(use)
  #here the original "ivar2index" function piece from propr package ends - this piece is intended to define a reference used for iqlr transformation (below)
  
  #here iqlr transformation begins
  message("Alert: Saving log-ratio transformed counts to @logratio.")
  logX <- log(counts)
  logSet <- logX[, use, drop = FALSE]
  ref <- rowMeans(logSet)
  lr <- sweep(logX, 1, ref, "-")
  
  return(lr)
}
################################




## upload and prepare data. more detailed comments below each line.

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





########### Creating log-transformed table of SMAGs counts

dataSMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs_for_SAGs_for_all.xlsx", sheet = 7, startRow = 1, colNames = TRUE) %>% data.table
dataSMAGs <- dataSMAGs[1:713,]

##########################RELATIVE ABUNDANCE CALCULATION OMITTED DUE TO THE IQLR TRANSFORMATION ADDED BELOW################
#calculating the relative abundance of SMAGs

#numeric_cols <- names(dataSMAGs)[sapply(dataSMAGs, is.numeric)]
#dataSMAGs_relabund <- dataSMAGs[, lapply(.SD, function(x) x/sum(x)), .SDcols = numeric_cols]
#dataSMAGs <- cbind(dataSMAGs[,1],dataSMAGs_relabund)
#calculating the relative abundance of SMAGs
##########################RELATIVE ABUNDANCE CALCULATION OMITTED DUE TO THE IQLR TRANSFORMATION ADDED BELOW################

dataSMAGs_t <- t(dataSMAGs)
colnames(dataSMAGs_t) <- dataSMAGs_t[1,]

dataSMAGs_t1 <-(as.data.frame(dataSMAGs_t))

dataSMAGs_t1 <- data.frame(dataSMAGs_t1[-1,])

#colnames(dataSMAGs_t1) <- NULL
#rownames(dataSMAGs_t1) <- NULL

dataSMAGs_t1 <- data.frame(apply(dataSMAGs_t1, 2, function(x) as.numeric(as.character(x))))


################################IQLR TRANSFORMATION FROM PROPR######################################

#dataSMAGs_t2 <- ivar2index(dataSMAGs_t1)

#OR

dataSMAGs_t2 <- propr(dataSMAGs_t1, # rows as samples, like it should be
                      metric = "rho", # or "phi", "phs", "cor", "vlr"
                      ivar = "iqlr", # or can use "iqlr" instead
                      #alpha = NA, # use to handle zeros
                      p = 100) # used by updateCutoffs

dataSMAGs_t2 <- dataSMAGs_t2@logratio

################################IQLR TRANSFORMATION FROM PROPR######################################

rownames(dataSMAGs_t2) <- colnames(dataSMAGs)[-1]

#dataSMAGs_t2$variable <- rownames(dataSMAGs_t2) I dont know why I made this line. Results in errors. probably I thought "transpose" function applied below requires varaible, but that was not true.
#dataSMAGs_t2 <- t(dataSMAGs_t2)

dataSMAGs_t2 <- as.data.frame(t(dataSMAGs_t2))

#dataSMAGs_t2 <- data.frame(apply(dataSMAGs_t2, 2, function(x) as.numeric(as.character(x))))

#rownames(dataSMAGs_t2) <- colnames(dataSMAGs)[-1]
dataSMAGs_t2$SMAG <- rownames(dataSMAGs_t2)

#########################################

########### Code used to create log-transformed table of V9 counts



#reading and reshaping V9 dataset
#dataV9 <- fread("18S_V9_1476_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.dereplicated.table (3)")

#dataV9 <- dataV9[,c(1, 5:1088)] # here ncol(dataV9) =1088; this line is to discardcolumns with total abundance and other miscellaneous abundance stats
##########################RELATIVE ABUNDANCE CALCULATION OMITTED DUE TO THE IQLR TRANSFORMATION ADDED BELOW################

#dataV9_t <- t(dataV9)
#colnames(dataV9_t) <- dataV9_t[1,]
#dataV9_t1 <- as.data.frame(dataV9_t)
#dataV9_t1 <- dataV9_t1[-1,]
#dataV9_t1 <- data.frame(apply(dataV9_t1, 2, function(x) as.numeric(as.character(x))))


################################IQLR TRANSFORMATION FROM PROPR######################################

#dataV9_t2 <- ivar2index(dataV9_t1)
#dataV9_t2_BACKUP <- dataV9_t2

################################IQLR TRANSFORMATION FROM PROPR######################################

#rownames(dataV9_t2) <- colnames(dataV9)[-1]
#dataV9_t2 <- as.data.frame(t(dataV9_t2))
#dataV9_t2$amplicon <-str_replace(rownames(dataV9_t2), "X", "")
#dataV9_t2_BACKUP <- dataV9_t2

#########################################
#dataV9_t21 <- as.data.table(dataV9_t2)

#dataV9_melted <- dataV9_t21[,.SD,.SDcols=c("amplicon",grep("TARA",colnames(dataV9_t21),value=T))] %>%
#  melt(id.vars="amplicon")
##########writing this to the PRICELESS FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#fwrite(dataV9_melted, "dataV9_meltedLOG-TRANSFORMED.tsv",sep="\t",row.names = T)
##########writing this to the PRICELESS FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

########### Code used to create log-transformed table of V9 counts; counts saved in "dataV9_meltedLOG-TRANSFORMED.tsv"; as it takes >1h to calculate the log-transformed table, we start directly from the file created before.

#dataV9_melted1 <- fread("dataV9_meltedLOG-TRANSFORMED.tsv") #--the old version. To avoid reading this dataset each time and save time



diffmetr <- c("phi","rho","phs","cor","vlr") #a vector that stores different metrics of proportionality measurement

####preparing constant variables from V9 dataset for the Evil Shuffle 

dataV9_melted1_prim <- fread("dataV9_meltedLOG-TRANSFORMED.tsv")

dat <-dat <-dcast(dataV9_melted1_prim, variable~amplicon)
rownames(dat) <- dat$variable

datx <- dat[,-1]


  ##############################################################################
  #++++++++++++++++++++++++XXXXXXXXXXXXXXXXXXX BEGIN THE EVIL SHUFFLE XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for (iteration in c(1:100000)) {
    print("BEGIN THE EVIL SHUFFLE iteration")
    print(iteration)
    print("!!!!!!!!")
    #iteration <- 1
    
    ###########################shuffling V9 data
    
    datxfreq <- randomizeMatrix(as.matrix(datx), null.model = c("frequency"), iterations = 1000) #==SHUFFLE STATIONS
    datxfreqrich <- randomizeMatrix(as.matrix(datxfreq), null.model = c("richness"), iterations = 1000)
    
    
    datshuffled <- as.data.frame(datxfreqrich)
    
    datshuffled$variable <- dat$variable
    
    
    dataV9_meltedx <- data.table(datshuffled) %>% melt(id.vars="variable")
    #just adding dummy column and fixing colnames and column positions to reconstitute dataV9_melted1 look
    dataV9_meltedx$V1 <- c(1:nrow(dataV9_meltedx))
    colnames(dataV9_meltedx) <- c("variable", "amplicon", "value", "V1")
    dataV9_meltedx <- dataV9_meltedx[, c(4, 2, 1, 3)]
    
    
    dataV9_melted1 <- dataV9_meltedx
    ###########################shuffling V9 data
    
    
    for (metr in 1:length(diffmetr)){
      metric <- diffmetr[metr]
      print(metric)
      
      
      
    
    
    
    for (d in c(3:length(colnames(vocabulary)))){
      #d <- 3
      print(colnames(vocabulary)[d])
      
      list_of_all <- c()
      sg <- colnames(vocabulary)[d]
      
      for (o in list_of_all_old ) {if(is.na(vocabulary[which(vocabulary$taxa_SMAGs == o), which(colnames(vocabulary)==colnames(vocabulary)[d])])!=TRUE){list_of_all <- c(list_of_all, o)}}
      
      list_of_all
      #uncomment to test the example barcode-MAG pairs
      
      #i <- "Attheya"
      
      #uncomment to test the example barcode-MAG pairs
      
      for(i in list_of_all[1:(length(list_of_all))]) {
        
        
        #subsetting SMAGs dataset by the taxogroup
        
        tokeep <- c(taxa_SMAGs_melted[which(taxa_SMAGs_melted$AA_Best_taxonomy_GENRE == i),2])
        
        #uncomment to test the example barcode-MAG pairs
        #tokeep <- c("TARA_ARC_108_MAG_00240","TARA_ARC_108_MAG_00247","TARA_ARC_108_MAG_00250","TARA_ARC_108_MAG_00255",
        #            "TARA_ARC_108_MAG_00257","TARA_ARC_108_MAG_00260","TARA_MED_95_MAG_00420","TARA_MED_95_MAG_00435",
        #           "TARA_PSE_93_MAG_00240","TARA_SOC_28_MAG_00054","TARA_SOC_28_MAG_00062")
        #uncomment to test the example barcode-MAG pairs
        
        #dataSMAGs_t2 <- dataSMAGs_t2[SMAG%in%tokeep] ##old solution
        
        dataSMAGs_t21 <- as.data.table(dataSMAGs_t2)[SMAG%in%tokeep]
        
        #subsetting SMAGs dataset by the taxogroup
        
        
        #reshaping the dataframe and matching station IDs with metagenome IDs
        
        dataSMAGs_melted <- melt(dataSMAGs_t21,id.vars="SMAG")
        
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
        
        #reshaping the dataframe and matching station IDs with metagenome IDs
        
        #writing out intermediate dataset of SMAGs abundances
        
        fwrite(dataSMAGs_df, paste0(metric, i, sg, "data_SMAGs_inv_ext_shuff.tsv"),sep="\t",row.names = T)
        
        
        #writing out intermediate dataset of SMAGs abundances
        
        ###reading V9 data
        
        dataV9 <- fread("18S_V9_1476_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.dereplicated.table (3)")
        dataV9 <- dataV9[,c(1, 5:1088)] 
        
        #dataV9_melted[,value2:=value/sum(value),by=variable]
        ####################################################################
        #the right newer version for the new SMAG subsets and for including SAGs
        taxogroup_i <- data.frame(vocabulary)[which(vocabulary == i)[1], d]
        amplicons <- datataxa[grep(taxogroup_i, datataxa$taxonomy),"amplicon"]
        
        #outdated fucking solution for non-extended data because of which I lost a day
        #taxogroup_i <- vocabulary[which(vocabulary == i)[1], colnames(vocabulary)[d]]
        #amplicons <- datataxa[grep(i, datataxa$taxonomy),"amplicon"]
        
        #not all amplicons present in datataxa are also present in dataV9. THis is exeptionally strange, and should be investigated in more detail later. I will treat it as artifact of analysis by now. but weird. the line below is a solution to treat such cases
        
        if (length(which(dataV9$amplicon %in% amplicons$amplicon))==0){next} else
          
          dataV9 <- dataV9[which(dataV9$amplicon %in% amplicons$amplicon),]
        
        #dataV9 <- dataV9[taxogroup==c(taxogroup_i)[1]&pid>90]
        
        x <- apply(dataV9[,.SD,.SDcols=grep("TARA",colnames(dataV9),value=T)],1,function(X) sum(X>0))
        dataV9 <- dataV9[x>100]
        
        if (length(dataV9$amplicon)==0){next} else
          
          dataV9_melted <- dataV9_melted1[amplicon%in%dataV9[,amplicon]]
        
        
        
        contextV9 <- fread("TARA_CONTEXT_95_SEQ_STATIONS_V9.tsv")
        dataV9_melted <- merge(dataV9_melted,contextV9[,list(Pangaea_sample_id,Sample.mat)],by.x="variable",by.y="Pangaea_sample_id")
        dataV9_melted[,variable:=NULL]
        dataV9_melted[,V1:=NULL]
        dataV9_df_relabund <- dcast(Sample.mat~amplicon,value.var="value",data=dataV9_melted,fun.aggregate = mean) %>%
          data.frame(row.names = "Sample.mat",check.names = F)
        
        
        
        
        fwrite(dataV9_df_relabund,paste0(metric, i,  sg, iteration, "_rand_","data_V9_relabund_inv_ext_shuff.tsv"),sep="\t",row.names = T)
        
        # writing out intermediate dataset of V9 abundances
        
        
        #######################Correlation estimation##############
        tmp1 <- dataSMAGs_df
        
        ##############PROPR PACKAGE EDITION STARTS HERE################
        tmp1 <- na.omit(tmp1)
        #deleting NAs; 0s will be automatically replaced by propr package function
        #tmp1[is.na(tmp1)==TRUE] <- 0
        tmp1 <- merge(as.data.frame(tmp1),as.data.frame(dataV9_df_relabund),by="row.names") %>%
          data.frame(row.names = "Row.names",check.names = F)
        #Stores relabund of all SMAGs and all barcodes(columns) in all samples(rows) 
        
        
        
        #####################################################################################################################################################3333
        ##############################proportionality calculation#####################################
        #metric <- "rho"
        lr <- as.matrix(tmp1)
        lrv <- propr:::lr2vlr(lr)
        
        if(metric == "rho"){
          mat <- propr:::lr2rho(lr)
        }else if(metric == "phi"){
          mat <- propr:::lr2phi(lr)
          #if(symmetrize) symRcpp(mat) # optionally force symmetry
        }else if(metric == "phs"){
          mat <- propr:::lr2phs(lr)
        }else if(metric == "cor"){
          mat <- stats::cor(lr)
        }else if(metric == "vlr"){
          mat <- lrv
        }else{
          stop("Provided 'metric' not recognized.")
        }
        pr<-mat
        
        
        #pr <- propr:::lr2rho(as.matrix(tmp1))
        colnames(pr) <- colnames(tmp1)
        rownames(pr) <- colnames(tmp1)
        
        labels <- propr:::labRcpp(ncol(as.matrix(tmp1)))
        result <-
          data.frame(
            "Partner" = colnames(pr)[labels[[1]]],#in original propr package code it was just "Partner" = labels[[1]] , had to be replaced
            "Pair" = colnames(pr)[labels[[2]]],#in original propr package code it was just "Pair" = labels[[2]] , had to be replaced
            "lrv" = propr:::lltRcpp(pr),
            "metric" = factor(metric),
            "alpha" = "Ignored for this time",
            "propr" = propr:::lltRcpp(pr)
          )
        
        ###small example to prove the replacement of "Partner" and "Pair" assignment is working as expected    
        # x <- data.frame(
        #  "A" = c (1,2,3),
        # "B" = c (1,2,3),
        #"C" = c (1,2,3))
        #  colnames(x) <- c("a","b","c")
        # rownames(x) <- c("a","b","c")
        #labels <- propr:::labRcpp(ncol(as.matrix(x)))
        #colnames(x)[labels[[1]]]
        
        
        #pr <- propr(tmp1, # rows as samples, like it should be
        #           metric = "rho", # or "phi", "phs", "cor", "vlr"
        #          ivar = "iqlr", # or can use "iqlr" instead
        #         #alpha = NA, # use to handle zeros
        #        p = 100) # used by updateCutoffs
        
        #res <- getResults(pr)
        res <- result
        #make a dataframe from the propr package output
        res <- res[grep("TARA|Metagenome_centric_SAG",res$Pair),]
        res <- res[grep("TARA|Metagenome_centric_SAG",res$Partner, invert = T),]
        #taking only barcode+MAG pairs
        
        #res[grep("TARA",res$Partner, invert = T),]
        #checking if it worked; uncomment if necessary
        
        res <- res[grep(".y",res$Partner, invert = T),]
        res <- res[grep(".x",res$Partner, invert = T),]
        # for some reason, function produces 3 outputs per each partner, adding .x , .y and nothing on the end of the partner name. Rho values do not differ between these three... Until I figure out the reason why it does this way, I delete them for now...
        
        
        #now merge
        data_cor <- merge(res,dataV9[,list(amplicon)],by.x="Partner",by.y="amplicon")
        data_cor <- merge(res,dataV9[,list(amplicon)],by.x="Partner",by.y="amplicon")
        if (length(levels(as.factor(data_cor$V2))) == 1){data_cor$V2<-dataV9[1,list(amplicon)]}
        data_cor <- data_cor[,c(1,2,6)]
        setnames(data_cor,c("md5sum","SMAG","rho"))
        #adds some info about the barcodes to the data_corr. ctotab is  total abundance of barcodes in this OTU.
        
        #data_cor[order(rho,decreasing = T)][1:10]
        #What was the reason for this line?.....
        
        fwrite(data_cor,paste0(metric, i, sg, iteration, "_rand_", "all_scores_relabund_inv_ext_shuff.tsv"),sep="\t")
        
        
        ########code used for common correlation follows in the comment below########
        #if(length(tmp1[grep("TARA",colnames(tmp1))]) == 1) next
        #data_cor <- cor(as.data.frame(tmp1[,grep("TARA",colnames(tmp1))]),as.data.frame(tmp1[,grep("TARA",colnames(tmp1),invert = #T)]),method="spearman",use="pairwise.complete.obs") %>%
        #data.table(keep.rownames = TRUE) %>%
        #melt(id.vars="rn")
        #############################################################################
        
        #produces correllation value matrix where for every barcode-MAG pair correlatio coefficient is stored(in rows)
        #data_cor <- merge(res,dataV9[,list(md5sum,pid,lineage,sequence,refs,ctotab)],by.x="Partner",by.y="md5sum")
        #setnames(data_cor,c("md5sum","SMAG","rho","metric", "alpha", "proper", "zeros", "pid","lineage","sequence","refs","abundance"))
        #adds some info about the barcodes to the data_corr. ctotab is  total abundance of barcodes in this OTU.
        
        #data_cor[order(rho,decreasing = T)][1:10]
        #What was the reason for this line?.....
        
        #fwrite(data_cor,paste0(i,"all_scores_relabund.tsv"),sep="\t")
        
        data_cor <- data.table(data_cor)
        #necessary after propr package function
        split(data_cor,data_cor[,"SMAG"]) %>%
          lapply(function(X){
            X[order(rho,decreasing = T),][1:3]
          }) %>% rbindlist() %>% fwrite(paste0(metric, i, sg, iteration, "_rand_", "best_3_scores_relabund_inv_ext_shuff.tsv"),sep="\t")
        #splits dataframe according to SMAG, and takes only top-3 correllation coefficients for every SMAG 
        
      }
    }
  }
}