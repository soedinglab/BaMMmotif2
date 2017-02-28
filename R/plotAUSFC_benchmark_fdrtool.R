#!/usr/bin/Rscript
#-----------------------------------------------------------------
#
# Documentation
#
#-----------------------------------------------------------------

# get the p-values for both positives and negatives from the
# .zoops.stats file;
# calculate false-discovery-rate (FDR) using modified fdrtool
# with a fixed eta0, which equals fold/(1+fold);
# plot the FDR vs. recall(sensitivity) curve;
# caluclate the area under the sensitivity-FDR curve (AUSFC)
# and plot the true-positive-rate(TPR) vs. false-positive-rate(FPR)
# curve and calculate the partial AUC from this curve;
# results are saved in a .bmscore file.

# examples for running this script:
# ./plotAUSFC_benchmark_fdrtool.R PATH_TO_zoops.stats_FILE BASENAME_OF_THE_FILE FACTOR
# ./plotAUSFC_benchmark_fdrtool.R /home/bamm_result/ JunD_motif_1 10

#-----------------------------
#
# R sources
#
#-----------------------------

# load "zoo" package for calculating AUPRC
library( zoo )

# fdrtool function from the same folder as this script
source( "fdrtool.R", chdir = TRUE )

#-----------------------------
#
# Parameter setting
#
#-----------------------------

args <- commandArgs(trailingOnly=TRUE)

# get the directory of the p-value files
dir <- args[1]
dataname <- args[2]
mfold <- as.numeric(args[3])

# read in p-values from file
stats <- read.table(paste(dir, dataname, ".zoops.stats", sep = "" ), skip=1 )
pvalues <- stats$V5

# avoid the rounding errors when p-value = 0 or p-value > 1
for(i in seq(1, length(pvalues))){
  if( pvalues[i] > 1 ){
    pvalues[i] = 1
  }
}

# estimate False Discovery Rates for Diverse Test Statistics
eta0set = mfold/(1+mfold)
png( file = paste( dir, dataname, '_FDRstat.png', sep = "" ) )
result = fdrtool( pvalues, statistic="pvalue",
                  plot=TRUE, color.figure=TRUE, verbose=TRUE, eta0set=eta0set )
dev.off()

# get the global fdr values and estimate of the weight eta0 of 
# the null component
fdr <- result$qval
eta0 <- result$param[3]

# calculate recall
len = length(fdr)
list <- seq(1, len )
recall <- ( 1 - fdr ) * list / ( 1 - eta0 ) / len

# modify the SF curve:
# reset recall to 1 when it is larger than 1
for(i in list){
  if( recall[i] > 1 ){
    for(rest in i:len){
      recall[rest] = 1
    }
    break
  }
}

# extend the SF curve till FDR reaches 0.5
recall <- append(recall, 1)
fdr <- append(fdr, 0.5)
range <- seq(1, len+1)

# limit the frame of the curve
left = 1
right = len+1
for(i in range){
  if( fdr[i] >= 0.0001 ){
    left = i
    recall[i] = 0
    break
  }
}

for(i in range){
  if(fdr[i] >= 0.5 ){
    right = i
    break
  }
}

range <- seq(left, right)

# plot fdr vs. recall curve
sum_area = log10(0.5) + 4
png( file = paste( dir, dataname, '_SFCurve.png', sep = "" ) )
plot(fdr[range], recall[range], log="x", 
     main=paste(dataname, "Sensitivity vs. FDR", sep="\n\n"),
     xlab="FDR", ylab="Sensitivity", xlim=c(0.0001,0.5), ylim=c(0,1),
     type='l', lwd=2.5,
     col="deepskyblue")
# compute  the area under the sensitivity-FDR curve(AUSFC):
if(right == left){
  ausfc = 0
} else {
  ausfc = sum(diff(log10(fdr[range]))*rollmean(recall[range],2)) / sum_area
}
# write the AUSFC on the plot with the precision of 4 digits:
text(0.0005, 0.95, paste( "AUSFC: ", round(ausfc, digits=4) ))
dev.off()

# plot TPR vs FPR
TP <- stats$V1
FP <- stats$V2
TPR <- TP / TP[length(TP)]
FPR <- FP / FP[length(FP)]
for(i in seq(1,length(FP))){
  if( FPR[i] >= 0.05 ){
    rbound = i
    break
  }
}
png( file = paste( dir, dataname, '_pTFCurve.png', sep = "" ) )
plot(FPR[1:rbound], TPR[1:rbound],
     main=paste(dataname, " FPR vs. TPR ", sep="\n\n"),
     xlab="FPR", ylab="TPR", xlim=c(0,0.05), ylim=c(0,1),
     type='l', lwd=2.5, col="deepskyblue")
# compute the area under the partical TP-FP curve(AUC5):
auc5 = sum(diff(FPR[1:rbound])*rollmean(TPR[1:rbound],2)) / 0.05

# write the AUC on the plot with the precision of 4 digits:
text(0.03, 0.1, paste( "pAUC: ", round(auc5, digits=4) ))
dev.off()


# write paramaters to a file
ofile = paste( dir, dataname, '.bmscore', sep = "" )
writeLines( c( "AUSFC: ", ausfc, 
               "pAUC: ", auc5), ofile )
