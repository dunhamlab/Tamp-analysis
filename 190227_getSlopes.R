#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## args
# args[1] = scripts directory
# args[2] = working directory
# args[3] = group (MN_R)

group <- unlist(strsplit(args[3],'_',fixed = TRUE))

library(reshape2)
library(grDevices)
source(paste(args[1],'/190227_CalculateSlopes.R',sep = ''))

# X is our data
X <- read.csv(paste(args[2],'/',args[3],'_log2.csv',sep=''),header=T,sep=",",row.names = NULL)
# X <- X[which(X$replicate=='B'),]
# X <- X[which(X$replicate=='A'),]

# get number of timepoints
cols <- ncol(X)
ngens <- (cols - 3)/2
gen <- melt(X[,c(1:3,(ngens+4):((2*ngens)+3))],id=c(1:3),variable.name="g",value.name = 'generation')
X <- melt(X[,c(1:3,4:(ngens+3))],id=c(1:3),variable.name="timepoint",value.name="growth")
# melted generations added to melted timepoints
X<-cbind(X,gen[,4:5])
# Reorder the data so that it is sorted by gene and X12mer
ord <- order(X$gene,X$mer,X$replicate)
X <- X[ord,] 

# #####################################################
# An example for you to try out the function for one gene
# #####################################################
# gene <- "YAL058W"
# slopes <- calculateSlope(X,gene,plot=TRUE)

# ######################################################
# If you want to do it for all the genes, you need to 
# create a for loop
# ######################################################
uniquegene <- unique(X$gene)
slopes <- NULL
se <- NULL
for(gene in uniquegene){
  temp <- calculateSlope(X,gene,dir=paste(args[2],'/',args[3],'_plots',sep=''),exp=args[3],plot=TRUE)
  temp2 <- summarizeSlopedensity(temp)
  # This contains the slope estimate and standard error
  # First element is the slope estimate
  # Second element is the standard error
  slopes <- c(slopes,temp2[1])
  se <- c(se, temp2[2])
}
slopes <- matrix(slopes,ncol=1)
se <- matrix(se,ncol=1)
results <- cbind(slopes,se)
colnames(results) <- c("slopes","standard error")
row.names(results) <- uniquegene
head(results)
write.csv(results, file = paste(args[2],'/',args[3],'_slopes.csv',sep=''))
