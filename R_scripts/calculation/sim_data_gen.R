library(stringr)
library(openxlsx)
library(data.table)
library(magrittr)
library(GGally)
library(propr)
library(picante)
library(plyr)
library(metaSPARSim)
require(GUniFrac)
require(vegan)
require(DESeq2)
##upload packages
library(stringr)
library(openxlsx)
library(data.table)
library(magrittr)
library(GGally)
library(propr)
##upload packages
##upload packages

##creaate modified function

GMPRX <- function(OTUmatrix, min_ct = 1, intersect_no = 1){
  #datatype check
  
  if(!(class(OTUmatrix)[1] %in% c("data.frame","matrix")))
    stop("Unknown datatype of object \"OTUmatrix\".")
  OTUmatrix  <- t(OTUmatrix)
  rownames(OTUmatrix)->SampleName
  if(ncol(OTUmatrix)<nrow(OTUmatrix))
    warning("Sample size is larger than OTU number. Check if samples are arranged in columns.")
  if(length(OTUmatrix[OTUmatrix<1 & OTUmatrix>0])>=0.1*ncol(OTUmatrix)*nrow(OTUmatrix))
    stop("More than 10% values are fractional, please check.")
  apply(OTUmatrix, MARGIN = 2, as.integer)->OTUmatrix
  
  
  GUniFrac:::gmpr(OTUmatrix,min_ct,intersect_no)->size.factor
  names(size.factor)<-SampleName
  size.factor[abs(size.factor-1)<1e-10]<-NA
  size.factor
}


##function for log-ratio transformation
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
##



################################
## set working directory

setwd("/home/dzavadska/simulation_pipeline/")


#Simulating V9 dataset 

#import sequencing depth vector; uncomment "#[1:5]" to run the mock example that doesnt consume so much memory 
#library sizes
varvect <- fread("V9_var_vector.tsv")$X#[1:5]

avgab <- fread("avg_abund_V9.tsv")$X#[1:5]
var_avgab <- fread("avg_abund_V9.tsv")$Y#[1:5]


#creating parameters

N_sample_cond_A <- length(varvect) 
param <- list()
param$intensity <- avgab #average abundances of V9 taxa
param$variability <- var_avgab #variance sensu the one in square
param$lib_size <- varvect #sequencing depth

for (i in c(76:80)) {
  #print("iteration", i)

  dataV9_melted <- NULL
#simulating dataset
metaSPARSim_param <- list(cond_A = param)
metaSPARSim_result <- metaSPARSim(metaSPARSim_param)
###
fwrite(metaSPARSim_result$count, paste0("/scratch/data1/dzavadska/sim_input","iteration_", i, "V9_big_simulation_priceless.tsv"))
###


##checking and verifying that metaSPARSim_result$counts is actually what we need
#str(metaSPARSim_result)
#dfSPARS <- as.data.frame(metaSPARSim_result$counts)
#sum(dfSPARS[,1])
#varvect[1]

#dfSPARSg <- as.data.frame(metaSPARSim_result$gamma)
#sum(dfSPARSg[,1])
#varvect[1]



##############recalculating SMAG simulation for raw counts

dataSMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs.xlsx", sheet = 7, startRow = 1, colNames = TRUE) %>% data.table

norm_total <- dataSMAGs[714,2:ncol(dataSMAGs)]

nbreadsSMAGs <- read.xlsx("Table_S04_distributions_nr_SMAGs.xlsx", sheet = 1, startRow = 1, colNames = TRUE) %>% data.table
nbreads <- nbreadsSMAGs[1:939,3]

dataSMAGs_raw <- data.frame(dataSMAGs)

total_raw <- as.numeric(norm_total)*as.numeric(nbreads$Nb.reads)/100

for (xyzab in 2:ncol(dataSMAGs[2:ncol(dataSMAGs)])) {
  dataSMAGs_raw[1:713,xyzab] <- dataSMAGs_raw[1:713,xyzab]*total_raw[xyzab-1]/as.numeric(norm_total)[xyzab-1] 
}



avgab_SMAGs <- apply(dataSMAGs_raw[,2:ncol(dataSMAGs)],1,function(X) mean(na.omit(as.numeric(X))))

varab_SMAGs <- apply(dataSMAGs_raw[,2:ncol(dataSMAGs)],1,function(X) var(na.omit(as.numeric(X))))


#this will also be sequencing depth

################################################
##try to estimate parameters from count table##
fwrite(dataSMAGs_raw[1:713,2:ncol(dataSMAGs_raw)], "SMAGs_raw_template.tsv")

#  normalizing data by GMPR, according to the instructions from MetaSPARSim package


tmp <- dataSMAGs_raw[1:713,2:ncol(dataSMAGs_raw)]
tmp <- tmp[,-which(tmp[1,]=="NaN")]

tmp <- tmp[,-which(total_raw == 0)]
tmp <- tmp[-which(avgab_SMAGs == 0),]


otu.tab <- t(tmp)
gmpr.size.factor <- GMPRX(otu.tab)


otu.tab.norm <- t(t(otu.tab)/gmpr.size.factor)
#str(otu.tab.norm)
tmp_norm <- t(otu.tab.norm)
#  normalizing data by GMPR, according to the instructions from MetaSPARSim package

indpops <- list(c(1:nrow(otu.tab)))
#row.names(tmp)<-dataSMAGs_raw[1:713,1]
params<-estimate_parameter_from_data(tmp, tmp_norm, indpops, perc_not_zeros = 0.01)

#now, the library size seems to large and the simulation runs for eternity. So, we divide the library size and intensity by 10, and variability by 100; also, library sizes in float format seem to take computation time, so we
# make lib sizes as integers.

modp <- list()
modp$intensity <- as.numeric(params[[1]][[1]])/10000
modp$variability <- as.numeric(params[[1]][[2]])/100000000
modp$lib_size <- as.integer(as.numeric(params[[1]][[3]])/10000)

metaSPARSim_modp <- list(cond_A = modp)

SMAGs_sim <- metaSPARSim(metaSPARSim_modp)
###
fwrite(SMAGs_sim$count, paste0("/scratch/data1/dzavadska/sim_input","iteration_", i, "MAG_SIMULATION_PRICELESS.tsv"))
################################################
#a dummy xmple

#dummy <- list()
#dummy$intensity <- as.numeric(params[[1]][[1]])[1:3]/1000
#dummy$variability <- as.numeric(params[[1]][[2]])[1:3]/1000000
#dummy$lib_size <- as.integer(as.numeric(params[[1]][[3]])[1:3]/1000)
#as.integer(as.numeric(params[[1]][[3]])[1:3]/10)
#metaSPARSim_dummy <- list(cond_A = dummy)
#SMAGs_dummy <- metaSPARSim(metaSPARSim_dummy)
#system.time(SMAGs_dummy <- metaSPARSim(metaSPARSim_dummy))



###############################################################################
###Distribution of intensities and adding spike-ins
#runs only on nisaba

simV9 <- fread(paste0("/scratch/data1/dzavadska/sim_input","iteration_", i, "V9_big_simulation_priceless.tsv"))
#str(simV9)

Y <- as.vector(t(simV9))


#create a spike-in vector
spikein <- sample (Y,size=1084)
spikein_ord <- spikein[order(as.numeric(spikein),decreasing = TRUE)]

######
#check if the number of 0s is approximately equal to the 0/non-0 ratio in the rest of the dataset
#frequency of 0s - ratio of 0 counts to non-0 counts
#null_to_nonnull <- length(which(Y!=0))/length(which(Y==0))
#null_to_nonnull = 0.3368987; length(Y) = 284789564; length(which(Y!=0)) = 71767017
#Y_nonnull <- Y[which(Y!=0)]
#summary(Y_nonull)
#levels(as.factor(Y_nonnull))

#length(which(spikein!=0))/length(which(spikein==0))
#null_to_nonnull 
#########


#creating a dataframe with spike-in
bind <- data.frame(matrix(NA, nrow = 1, ncol = 1084))
bind[1,] <- spikein_ord

simV9_spike1 <- rbind(simV9, bind, use.names=FALSE)
#creating a dataframe with spike-in

fwrite(simV9_spike1, paste0("/scratch/data1/dzavadska/sim_input","iteration_", i, "simV9_spike1_test,tsv"), sep = "\t")



#finding maximum values for each row 
#checking that the maxima for each column are by the magnitude of order larger comared to the spike-in values ==> no real need for substracting anything (however it can be done if badly requested)
#apply(simV9_spike1[1:262721,1:30],2,function(X) max(X))

simSMAG <- fread(paste0("/scratch/data1/dzavadska/sim_input","iteration_", i, "MAG_SIMULATION_PRICELESS.tsv"))

sample_V9_sum <- as.numeric(apply(simV9_spike1,2,function(X) sum(X)))
sample_SMAG_sum <- as.numeric(apply(simSMAG,2,function(X) sum(X)))

SMAG_spikein <- (as.numeric(simV9_spike1[262722,1:ncol(simSMAG)])*sample_SMAG_sum)/sample_V9_sum[1:ncol(simSMAG)]


#creating a dataframe with spike-in
bind_SMAG <- data.frame(matrix(NA, nrow = 1, ncol = ncol(simSMAG)))
bind_SMAG[1,] <- SMAG_spikein

simSMAG_spike1 <- rbind(simSMAG, bind_SMAG, use.names=FALSE)
#creating a dataframe with spike-in
fwrite(simSMAG_spike1, paste0("/scratch/data1/dzavadska/sim_input","positive","iteration_", i, "simSMAG_spike1_test,tsv"), sep = "\t")

#negative spike-in

neg <- as.numeric(simV9_spike1[262722,1:ncol(simSMAG)])
negative <- neg[order(as.numeric(neg),decreasing = FALSE)]

SMAG_spikein_n <- (negative*sample_SMAG_sum)/sample_V9_sum[1:ncol(simSMAG)]


#creating a dataframe with spike-in
bind_SMAG_n <- data.frame(matrix(NA, nrow = 1, ncol = ncol(simSMAG)))
bind_SMAG_n[1,] <- SMAG_spikein_n

simSMAG_spike1_n <- rbind(simSMAG, bind_SMAG_n, use.names=FALSE)
#creating a dataframe with spike-in
fwrite(simSMAG_spike1_n, paste0("/scratch/data1/dzavadska/sim_input","negative","iteration_", i, "simSMAG_spike1_test,tsv"), sep = "\t")



######
#check if the number of 0s is approximately equal to the 0/non-0 ratio in the rest of the dataset
#frequency of 0s - ratio of 0 counts to non-0 counts
#X <- as.vector(t(simSMAG))
#null_to_nonnullX <- length(which(X!=0))/length(which(X==0))
#null_to_nonnullX = 5.949712; length(Y) = 647463; length(which(X!=0)) = 554299
#X_nonnull <- X[which(X!=0)]
#summary(Y_nonull)
#levels(as.factor(Y_nonnull))

#length(which(spikein!=0))/length(which(spikein==0))
#null_to_nonnullX 
#########



###############################################
#LOG-RATIO TRANSFORMATION OF V9 dataset
###############################################




#reading and reshaping V9 dataset
dataV9 <- fread(paste0("/scratch/data1/dzavadska/sim_input", "iteration_", i, "simV9_spike1_test,tsv"))
dataV9$amplicon <- paste0("Y",rownames(dataV9))
dataV9 <- dataV9[,c(1085,1:1084)]
##########################RELATIVE ABUNDANCE CALCULATION OMITTED DUE TO THE IQLR TRANSFORMATION ADDED BELOW################

dataV9_t <- t(dataV9)
colnames(dataV9_t) <- dataV9_t[1,]
dataV9_t1 <- as.data.frame(dataV9_t)
dataV9_t1 <- dataV9_t1[-1,]
dataV9_t1 <- data.frame(apply(dataV9_t1, 2, function(x) as.numeric(as.character(x))))


################################IQLR TRANSFORMATION FROM PROPR######################################

dataV9_t2 <- ivar2index(dataV9_t1)
dataV9_t2_BACKUP <- dataV9_t2
fwrite(dataV9_t2_BACKUP, paste0("/scratch/data1/dzavadska/sim_input", "iteration_", i, "log-transfV9_simulation.tsv"),sep="\t",row.names = T)
################################IQLR TRANSFORMATION FROM PROPR######################################

rownames(dataV9_t2) <- colnames(dataV9)[-1]
dataV9_t2 <- as.data.frame(t(dataV9_t2))
dataV9_t2$amplicon <-str_replace(rownames(dataV9_t2), "X", "")
dataV9_t2_BACKUP <- dataV9_t2
#fwrite(dataV9_t2_BACKUP, "log-transfV9_simulation_backup2.tsv",sep="\t",row.names = T)
#########################################
dataV9_t21 <- as.data.table(dataV9_t2)

dataV9_melted <- dataV9_t21[,.SD,.SDcols=c("amplicon",grep("cond",colnames(dataV9_t21),value=T))] %>%  melt(id.vars="amplicon")
##########writing this to the PRICELESS FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
fwrite(dataV9_melted, paste0("/scratch/data1/dzavadska/sim_input", "iteration_", i, "dataV9_meltedLOG-TRANSFORMED_simulation.tsv"),sep="\t",row.names = T)
##########writing this to the PRICELESS FILE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


}









