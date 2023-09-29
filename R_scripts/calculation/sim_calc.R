## WARNING. the original script was run on nder Linux 5.15.0-46-generic, Ubuntu 22.04.1 with x86-64 architecture, on a machine with 32 CPU and 125Gb RAM, with R version 4.1.2 (2021-11-01). This is a computationally demanding which is likely not supposed to run on the local computers.


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

setwd("/scratch/data1/dzavadska/sim_out")

## set working directory



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


vocabulary <- fread("/scratch/data1/dzavadska/sim_out/taxa_counts.tsv",sep="\t") %>% data.frame()

list_of_all <- vocabulary[,1]
diffmetr <- c("phi","rho","phs","cor","vlr")

#######################FOR POSITIVE CORRESPONDANCE (ONLY, BY NOW)#####################

########loop for iteration
for(iteration in c(63:100)) {
print(iteration)
########loop for iteration
  #reading input data (different for each iteration)
  #reading the simulation data, and adjusting the name columns
  dataV9_pr <- fread(paste0("/scratch/data1/dzavadska/sim_input","iteration_", iteration, "simV9_spike1_test,tsv"))
  dataSMAGs <- fread(paste0("/scratch/data1/dzavadska/sim_input","positive","iteration_", iteration, "simSMAG_spike1_test,tsv")) %>% data.table
  
  dataSMAGs_t <- t(dataSMAGs)
  colnames(dataSMAGs_t) <- dataSMAGs_t[1,]
  
  dataSMAGs_t1 <-(as.data.frame(dataSMAGs_t))
  
  dataSMAGs_t1 <- data.frame(dataSMAGs_t1[-1,])
  
  dataSMAGs_t1 <- data.frame(apply(dataSMAGs_t1, 2, function(x) as.numeric(as.character(x))))
  ################################IQLR TRANSFORMATION FROM PROPR######################################
  dataSMAGs_t2 <- propr(dataSMAGs_t1, # rows as samples, like it should be
                        metric = "rho", # or "phi", "phs", "cor", "vlr"
                        ivar = "iqlr", # or can use "iqlr" instead
                        #alpha = NA, # use to handle zeros
                        p = 100) # used by updateCutoffs
  
  dataSMAGs_t2 <- dataSMAGs_t2@logratio
  
  ################################IQLR TRANSFORMATION FROM PROPR######################################
  rownames(dataSMAGs_t2) <- colnames(dataSMAGs)[-1]
  
  dataSMAGs_t2 <- as.data.frame(t(dataSMAGs_t2))
  
  dataSMAGs_t2$SMAG <- rownames(dataSMAGs_t2)
  #########################################
  
  
  
  #dataV9_melted1 <- fread(paste0("/scratch/data1/dzavadska/sim_input", "iteration_", iteration, "dataV9_meltedLOG-TRANSFORMED_simulation.tsv"))
  
  
  
#NOW START ACTUALLY LOOPING WITHIN THE ITERATION  
  
  ########loop for SMAGs amount
  for(i in list_of_all[1:(length(list_of_all))]) {
  ########loop for SMAGs amount
  print(i)
  pseudo_i <- which(vocabulary[,1]==i)
  num_of_SMAGs <-  vocabulary[pseudo_i,2]
  tokeep <- c(sample(1:(nrow(dataSMAGs_t2)-1),size= num_of_SMAGs-1))
  
  
  #creating a big subset of V9s at the different taxlevels; smaller subsets to be created from the bigger subsets from higher taxlevel
 
  #print(nrow(dataV9_pr))
  
  #note also that in a single case of Geminigera there is an issue when after subssetting V9 for supergroup they are less in the supergroup than in the taxogroup 1. So an exception is witten for this occasion. Normally, this in theory can happen and it is OK, on practice it happened only once in this particular dataset.
  
  if(i == "Geminigera"){
    V9_supergroup <- c(sample (1:(nrow(dataV9_pr)-1),size=vocabulary[pseudo_i,3]-1),  nrow(dataV9_pr))
    stupid_geminigera <- c(1:(nrow(dataV9_pr)-1))[which(1:(nrow(dataV9_pr)-1) %in% V9_supergroup == FALSE)]
    
    if (vocabulary[pseudo_i,4]-1==-1) {V9_tax_1 <- integer(0)} else { V9_tax_1 <- c(sample (stupid_geminigera,size=vocabulary[pseudo_i,4]-vocabulary[pseudo_i,3]),V9_supergroup) }
     V9_tax_2 <-V9_tax_1
     
    if (vocabulary[pseudo_i,6]-1==-1) {V9_tax_middle <- integer(0)} else { V9_tax_middle <- c(sample (V9_tax_2,size=vocabulary[pseudo_i,6]-1),  nrow(dataV9_pr)) }
    if (vocabulary[pseudo_i,7]-1==-1) {V9_genus <- integer(0)} else { V9_genus <- c(sample (V9_tax_middle,size=vocabulary[pseudo_i,7]-1),  nrow(dataV9_pr))}
  } else{
    
     V9_supergroup <- c(sample (1:(nrow(dataV9_pr)-1),size=vocabulary[pseudo_i,3]-1),  nrow(dataV9_pr))
  if (vocabulary[pseudo_i,4]-1==-1) {V9_tax_1 <- integer(0)} else { V9_tax_1 <- c(sample (V9_supergroup,size=vocabulary[pseudo_i,4]-1),  nrow(dataV9_pr)) }
  if (vocabulary[pseudo_i,5]-1==-1) {V9_tax_2 <- integer(0)} else { V9_tax_2 <- c(sample (V9_tax_1,size=vocabulary[pseudo_i,5]-1),  nrow(dataV9_pr)) }
  if (vocabulary[pseudo_i,6]-1==-1) {V9_tax_middle <- integer(0)} else { V9_tax_middle <- c(sample (V9_tax_2,size=vocabulary[pseudo_i,6]-1),  nrow(dataV9_pr)) }
  if (vocabulary[pseudo_i,7]-1==-1) {V9_genus <- integer(0)} else { V9_genus <- c(sample (V9_tax_middle,size=vocabulary[pseudo_i,7]-1),  nrow(dataV9_pr))}
  
    
  }
 
 
  
  ########loop for taxogroup level
  for (d in c(3:length(colnames(vocabulary)))){
  ########loop for taxogroup level
    if (d == 3) {amplicons <- V9_supergroup} else if (d == 4) {amplicons <- V9_tax_1}  else if (d == 5) {amplicons <- V9_tax_2}  else if (d == 6) {amplicons <- V9_tax_middle}  else if (d == 7) {amplicons <- V9_genus}
    
    if (length(amplicons) == 0) {next}
    
    print(colnames(vocabulary)[d])
    sg <- colnames(vocabulary)[d]
    
    
    ########loop for metrics level    
      for (metr in 1:length(diffmetr)){
    ########loop for metrics level 
        metric <- diffmetr[metr]
       # print(metric)
        
        
        
        #subsetting SMAGs dataset by the taxogroup -- IN THIS CASE, TAXOGROUP IS SUBSTITUTED BY THE NUMBER (MIN, MAX AND AVERAGE) OF SMAGS IN A GENUS
        #tokeep <- c(sample (1:(nrow(dataSMAGs_t2)-1),size= i-1))
        
        dataSMAGs_t21 <- as.data.table(dataSMAGs_t2)[c(tokeep, nrow(dataSMAGs_t2)),]
        dataSMAGs_t21$SMAG <- paste0("TARA",rownames(dataSMAGs_t21)) #creating a "mock variable" that will mimic the "SMAG" column normally containing SMAG names
        ##NOTE:!!!!!!!!!!!!!the "spike-in" SMAG is always the last row!!!!!!!!!!!!!! 
        #subsetting SMAGs dataset by the taxogroup-- IN THIS CASE, TAXOGROUP IS SUBSTITUTED BY THE NUMBER (MIN, MAX AND AVERAGE) OF SMAGS IN A GENUS
        
        #reshaping the dataframe and matching station IDs with metagenome IDs(piece of matching station IDs with metagenome IDs omitted because there is no need to introduce it for simulation -- station IDs already match)
        dataSMAGs_melted <- melt(dataSMAGs_t21,id.vars="SMAG")
        
        dataSMAGs_df <- dcast(variable~SMAG,value.var="value",data=dataSMAGs_melted,fun.aggregate = mean) %>%
          data.frame(row.names = "variable")
        #reshaping the dataframe (piece of matching station IDs with metagenome IDs omitted because there is no need to introduce it for simulation -- station IDs already match)
        
        #writing out intermediate dataset of SMAGs abundances
        fwrite(dataSMAGs_df, paste0("/scratch/data1/dzavadska/sim_out/",metric, i, sg, iteration, "data_SMAGs_simulation_matr.tsv"),sep="\t",row.names = T)
        #writing out intermediate dataset of SMAGs abundances
        
        ###reading V9 data
        #dataV9 <- fread("simV9_spike1_test,tsv") --- moved outside of the loop to save calculation time
        dataV9_pr$amplicon <- paste0("Y",rownames(dataV9_pr))
        dataV9 <- dataV9_pr[,c(1085,1:1084)]
        
        ####################################################################
        
        dataV9 <- dataV9[amplicons,]
        
        x <- apply(dataV9[,.SD,.SDcols=grep("cond",colnames(dataV9),value=T)],1,function(X) sum(X>0))
        dataV9 <- dataV9[x>100]
        
        dataV9_melted1 <- fread(paste0("/scratch/data1/dzavadska/sim_input", "iteration_", iteration, "dataV9_meltedLOG-TRANSFORMED_simulation.tsv"))
        
        dataV9_melted <- dataV9_melted1[amplicon%in%dataV9[,amplicon]]
        
        #####! ! ! ! ! ! ! ! ! ! ! ! !
        #removing variable of a big size!!!
        rm(dataV9_melted1)
        #removing variable of a big size!!!
        #####! ! ! ! ! ! ! ! ! ! ! ! !
        
        dataV9_melted[,V1:=NULL]
        dataV9_df_relabund <- dcast(variable~amplicon,value.var="value",data=dataV9_melted,fun.aggregate = mean) %>%
          data.frame(row.names = "variable",check.names = F)
        
        fwrite(dataV9_df_relabund,paste0("/scratch/data1/dzavadska/sim_out/",metric, i,  sg, iteration, "data_V9_relabund_simulation_matr.tsv"),sep="\t",row.names = T)
        
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
        res <- res[grep("TARA",res$Pair),]
        res <- res[grep("TARA",res$Partner, invert = T),]
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
        
        fwrite(data_cor,paste0("/scratch/data1/dzavadska/sim_out/",metric, i, sg,  iteration,"all_scores_relabund_simulation_matr.tsv"),sep="\t")
        
        
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
          }) %>% rbindlist() %>% fwrite(paste0("/scratch/data1/dzavadska/sim_out/",metric, i, sg, iteration, "best_3_scores_relabund_simulation_matr.tsv"),sep="\t")
        #splits dataframe according to SMAG, and takes only top-3 correllation coefficients for every SMAG 
        
      }}}

  
  }
