#!/usr/bin/Rscript
#-----------------------------------------------------------------
#
# Documentation
#
#-----------------------------------------------------------------

# calculate FDR, precision, recall and occurrences for both MOPS 
# and ZOOPS models by taking p-values of log odds scores froom 
# positive sequences and using "fdrtool" library

#-----------------------------
#
# R sources
#
#-----------------------------

# load "fdrtool" library
if (!require("fdrtool")){
  install.packages("fdrtool", repos="http://cran.rstudio.com/")
}

# load "zoo" package for calculating AUPRC
if (!require("zoo")){
  install.packages("zoo", repos="http://cran.rstudio.com/")
}

#-----------------------------
#
# Parameter setting
#
#-----------------------------

args <- commandArgs(trailingOnly=TRUE)

# get the directory to p-value files
dir <- args[1]
dataname <- args[2]

# local tests:
#dir = c("/home/wanwan/benchmark/JunD/EM_FDR_new/")
#dataname = c("JunD_motif_1")

# read p-values from the files
pvalues_zoops <- unlist( read.table(paste(dir, dataname, ".zoops.pvalues", sep = "" )))
pvalues_mops <- unlist( read.table(paste(dir, dataname, ".mops.pvalues", sep = "" )))

# estimate False Discovery Rates for Diverse Test Statistics
png( file = paste( dir, dataname, '_FDRstat_ZOOPS.png', sep = "" ) )
result_zoops = fdrtool( pvalues_zoops, statistic="pvalue",
        plot=TRUE, color.figure=TRUE, verbose=TRUE )
dev.off()

png( file = paste( dir, dataname, '_FDRstat_MOPS.png', sep = "" ) )
result_mops = fdrtool( pvalues_mops, statistic="pvalue",
        plot=TRUE, color.figure=TRUE, verbose=TRUE )
dev.off()

# get the global fdr values and estimate of the weight eta0 of 
# the null component
fdr_zoops <- result_zoops$qval
eta0_zoops <- result_zoops$param[3]

fdr_mops <- result_mops$qval
eta0_mops <- result_mops$param[3]

# calculate precision and recall
len_zoops = length(fdr_zoops)
list_zoops <- seq(1, len_zoops )
precision_zoops <- 1 - fdr_zoops
recall_zoops <- ( 1 - fdr_zoops ) * list_zoops / ( 1 - eta0_zoops ) / len_zoops

len_mops = length(fdr_mops)
list_mops <- seq(1, len_mops )
precision_mops <- 1 - fdr_mops
recall_mops <- ( 1 - fdr_mops ) * list_mops / ( 1 - eta0_mops ) / len_mops

# calculate the occurrence fraction and occurrences per sequence
occ_frac = 1 - eta0_zoops
occ_mult = (1 - eta0_mops) * len_mops / len_zoops

# # plot 1- fdr vs. l for testing
# png( file = paste( dir, dataname, '_TPvsL.png', sep = "" ) )
# plot( list_mops, precision_mops * list_mops,
#      main=paste(dataname,"TP vs. l for MOPS model",sep="\n\n"), 
#      xlab="l", ylab="(1-Fdr)*l",
#      type='l', lwd=2.5,
#      col="deepskyblue",
#      cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
# dev.off()

# plot the precision-recall curves
# for ZOOPS:
png( file = paste( dir, dataname, '_PRcurve_ZOOPS.png', sep = "" ) )
plot(recall_zoops, precision_zoops,
     main=paste(dataname,"Precision vs. Recall curve for ZOOPS model",sep="\n\n"), 
     xlab="recall", ylab="precision", xlim=c(0,1), ylim=c(0,1),
     type='l', lwd=2.5,
     col="deepskyblue")
# compute  the area under the precision-recall curve(AUPRC):
auprc_zoops = sum(diff(recall_zoops[list_zoops])*rollmean(precision_zoops[list_zoops],2))
# write AUPRC on the plot:
text(0.8, 0.5, paste( "AUPRC = ", round(auprc_zoops, digits=3) ))
dev.off()

# for MOPS:
png( file = paste( dir, dataname, '_PRcurve_MOPS.png', sep = "" ) )
plot(recall_mops, precision_mops,
     main=paste(dataname,"Precision vs. Recall curve for MOPS model",sep="\n\n"), 
     xlab="recall", ylab="precision", xlim=c(0,1), ylim=c(0,1),
     type='l', lwd=2.5,
     col="deepskyblue")
# compute  the area under the precision-recall curve(AUPRC):
auprc_mops = sum(diff(recall_mops[list_mops])*rollmean(precision_mops[list_mops],2))
# write AUPRC on the plot:
text(0.8, 0.5, paste( "AUPRC = ", round(auprc_mops, digits=3) ))
dev.off()
     

# plot the log odd scores for checking
# read in tables and skip the first line
logodds_zoops <- read.table(paste(dir, dataname, ".zoops.logOdds", sep = "" ), skip=1)
logodds_mops <- read.table(paste(dir, dataname, ".mops.logOdds", sep = "" ), skip=1) 

png( file = paste( dir, dataname, '_LogOdds_ZOOPS.png', sep = "" ) )
hist(logodds_zoops$V1,
     main=paste(dataname,"log odds score distribution for ZOOPS model",sep="\n\n"), 
     xlab="Log odds scores", ylab="Density",
     col=rgb(1,0,0,1/4), density=NULL, breaks=50)
hist(logodds_zoops$V2, add=T, col=rgb(0,0,0,1/4), density=NULL, breaks=50)
dev.off()

png( file = paste( dir, dataname, '_LogOdds_MOPS.png', sep = "" ) )
hist(logodds_mops$V1,
     main=paste(dataname,"log odds score distribution for MOPS model",sep="\n\n"), 
     xlab="Log odds scores", ylab="Density",
     col=rgb(1,0,0,1/4), density=NULL, breaks=50)
hist(logodds_mops$V2, add=T, col=rgb(0,0,0,1/4), density=NULL, breaks=50)
dev.off()

# plot the p-values for checking
png( file = paste( dir, dataname, '_Pvalues_ZOOPS.png', sep = "" ) )
hist(pvalues_zoops,
     main=paste(dataname,"p-value distribution for ZOOPS model",sep="\n\n"), 
     xlab="p-values", ylab="Density",
     col="skyblue", breaks=20, density = NULL,
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
dev.off()

png( file = paste( dir, dataname, '_Pvalues_MOPS.png', sep = "" ) )
hist(pvalues_mops,
     main=paste(dataname,"p-value distribution for MOPS model",sep="\n\n"), 
     xlab="p-values", ylab="Density",
     col="skyblue", breaks=20, density = NULL,
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5)
dev.off()

# only for testing:
# plot the recall-FDR curves from calculatePR() function of the c++ code
# read in tables from files, and skip the first line
zoops_stat <- read.table(paste(dir, dataname, ".zoops.stats", sep = "" ), skip=1)
l_zoops_stat <- seq(1, length(zoops_stat$V1) )
fdr_cpp_zoops <- zoops_stat$V3
for(i in l_zoops_stat){
  if( fdr_cpp_zoops[i] < 0.000001 ){
    fdr_cpp_zoops[i] = 0.000001
  }
}
recall_cpp_zoops <- zoops_stat$V4
png( file = paste( dir, dataname, '_SFcurve_ZOOPS_cpp.png', sep = "" ) )
plot( log10(fdr_cpp_zoops), recall_cpp_zoops, 
     main=paste(dataname,"Recall vs. FDR curve for ZOOPS model(c++)",sep="\n\n"), 
     xlab="log10(FDR)", ylab="Recall", xlim=c(-4,0), ylim=c(0,1),
     type='l', lwd=2.5, 
     col="deepskyblue")
# calculate the area under the SF curve:
ausfc_zoops = sum(diff(log10(fdr_cpp_zoops)[l_zoops_stat])*rollmean(recall_cpp_zoops[l_zoops_stat],2))
# print the AUSFC on the plot with the precision of 3-digits:
text(-1, 0.5, paste( "AUSFC = ", round(ausfc_zoops, digits=3) ))
dev.off()

mops_stat <- read.table(paste(dir, dataname, ".mops.stats", sep = "" ), skip=1)
l_mops_stat <- seq(1, length(mops_stat$V1) )
fdr_cpp_mops <- mops_stat$V3
for(i in l_mops_stat){
  if( fdr_cpp_mops[i] < 0.000001 ){
    fdr_cpp_mops[i] = 0.000001
  }
}
recall_cpp_mops <- mops_stat$V4
png( file = paste( dir, dataname, '_SFcurve_MOPS_cpp.png', sep = "" ) )
plot(log10(fdr_cpp_mops), recall_cpp_mops,
     main=paste(dataname,"Recall vs. FDR curve for MOPS model(c++)",sep="\n\n"), 
     xlab="log10(FDR)", ylab="Recall", xlim=c(-4,0), ylim=c(0,1),
     type='l', lwd=2.5,
     col="deepskyblue")
# calculate the area under the SF curve:
ausfc_mops = sum(diff(log10(fdr_cpp_mops)[l_mops_stat])*rollmean(recall_cpp_mops[l_mops_stat],2))
# print the AUSFC on the plot with the precision of 3-digits:
text(-1, 0.5, paste( "AUSFC = ", round(ausfc_mops, digits=3) ))
dev.off()

# write paramaters to a file
ofile = paste( dir, dataname, '.FDRtool', sep = "" )
writeLines(c("AUPRC(MOPS): ", auprc_mops, 
             "AUPRC(ZOOPS): ", auprc_zoops, 
             "fraction occurrence: ", occ_frac, 
             "multiple occurrence: ", occ_mult,
             "eta0(MOPS): ", eta0_mops,
             "eta0(ZOOPS): ", eta0_zoops), ofile)

