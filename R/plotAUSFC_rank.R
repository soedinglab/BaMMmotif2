#!/usr/bin/Rscript
#-----------------------------------------------------------------
#
# Documentation
#
#-----------------------------------------------------------------

# calculate FDR, precision, recall and occurrences for both MOPS 
# and ZOOPS models by taking p-values of log odds scores from 
# positive sequences and using "fdrtool" library;
# plot FDR vs. recall(sensitivity) curve;
# caluclate the area under the sensitivity-FDR curve (AUSFC);
# results are saved in a .rankscore file;

# examples for running this script:
# ./plotAUSFC_benchmark_fdrtool.R PATH_TO_zoops.stats_FILE BASENAME_OF_THE_FILE
# ./plotAUSFC_benchmark_fdrtool.R /home/bamm_result/ JunD_motif_1

#-----------------------------
#
# R sources
#
#-----------------------------

# load "fdrtool" library for calculating FDR and recall
library( fdrtool )

# load "zoo" package for calculating AUPRC
library( zoo )

#-----------------------------
#
# Parameter setting
#
#-----------------------------

args <- commandArgs(trailingOnly=TRUE)

# get the directory of the p-value files
dir <- args[1]
dataname <- args[2]

# read in p-values from the files
pvalues_zoops <- unlist( read.table(paste(dir, dataname, ".zoops.pvalues", sep = "" )))
pvalues_mops <- unlist( read.table(paste(dir, dataname, ".mops.pvalues", sep = "" )))

# avoid the rounding errors when p-value = 0 or p-value > 1
for(i in seq(1, length(pvalues_zoops))){
  if( pvalues_zoops[i] > 1 ){
    pvalues_zoops[i] = 1
  }
}
for(i in seq(1, length(pvalues_mops))){
  if( pvalues_mops[i] > 1 ){
    pvalues_mops[i] = 1
  }
}

# estimate False Discovery Rates for Diverse Test Statistics
result_zoops = fdrtool( pvalues_zoops, statistic="pvalue",
                        plot=FALSE, color.figure=TRUE, verbose=TRUE )

result_mops = fdrtool( pvalues_mops, statistic="pvalue",
                       plot=FALSE, color.figure=TRUE, verbose=TRUE )

# get the global fdr values and estimate of the weight eta0 of 
# the null component
fdr_zoops <- result_zoops$qval
eta0_zoops <- result_zoops$param[3]
fdr_mops <- result_mops$qval
eta0_mops <- result_mops$param[3]

# avoid eta0 = 1
if(eta0_zoops==1){
  eta0_zoops=0.9999999
}
if(eta0_mops==1){
  eta0_mops=0.9999999
}

# calculate recall
len_zoops = length(fdr_zoops)
list_zoops <- seq(1, len_zoops )
recall_zoops <- ( 1 - fdr_zoops ) * list_zoops / ( 1 - eta0_zoops ) / len_zoops

len_mops = length(fdr_mops)
list_mops <- seq(1, len_mops )
recall_mops <- ( 1 - fdr_mops ) * list_mops / ( 1 - eta0_mops ) / len_mops

# modify the SF curve:
# reset recall to 1 when it is larger than 1
for(i in list_mops){
  if( recall_mops[i] > 1 ){
    for(rest in i:len_mops){
      recall_mops[rest] = 1
    }
    break
  }
}

for(i in list_zoops){
  if( recall_zoops[i] > 1 ){
    for(rest in i:len_zoops){
      recall_zoops[rest] = 1
    }
    break
  }
}

# extend the SF curve till FDR reaches 0.5
recall_zoops <- append(recall_zoops, 1)
fdr_zoops <- append(fdr_zoops, 0.5)
recall_mops <- append(recall_mops, 1)
fdr_mops <- append(fdr_mops, 0.5)

list_zoops <- seq(1, len_zoops+1)
list_mops <- seq(1, len_mops+1)

# limit the frame of the curve
left_zoops = 1
right_zoops = len_zoops+1
left_mops = 1
right_mops = len_mops+1
for(i in list_zoops){
  if( fdr_zoops[i] >= 0.0001 ){
    left_zoops = i
    break
  }
}

for(i in list_zoops){
  if(fdr_zoops[i] >= 0.5 ){
    right_zoops = i
    break
  }
}

list_zoops <- seq(left_zoops, right_zoops)

for(i in list_mops){
  if(fdr_mops[i] >= 0.0001 ){
    left_mops = i
    break
  }
}
for(i in list_mops){
  if(fdr_mops[i] >= 0.5 ){
    right_mops = i
    break
  }
}

list_mops <- seq(left_mops, right_mops)

sum_area = log10(0.5) + 4

# calculate the occurrence fraction and occurrences per sequence
occ_frac = 1 - eta0_zoops
occ_mult = (1 - eta0_mops) * len_mops / len_zoops

# plot fdr vs. recall curve
png( file = paste( dir, dataname, '_SFCurve_ZOOPS.png', sep = "" ) )
plot(fdr_zoops[list_zoops], recall_zoops[list_zoops], log="x", 
      main=paste(dataname, "Sensitivity vs. FDR for ZOOPS model", sep="\n\n"),
      xlab="FDR", ylab="Sensitivity", xlim=c(0.0001,0.5), ylim=c(0,1),
      type='l', lwd=2.5,
      col="deepskyblue")
# compute  the area under the sensitivity-FDR curve(AUSFC):
if(right_zoops == left_zoops){
  ausfc_zoops = 0
} else {
  ausfc_zoops = sum(diff(log10(fdr_zoops[list_zoops]))*rollmean(recall_zoops[list_zoops],2)) / sum_area
}

# write the AUSFC on the plot with the precision of 3 digits:
text(0.0005, 0.95, paste( "AUSFC: ", round(ausfc_zoops, digits=4) ))
# write the fraction of occurrence on the plot with the precision of 3 digits:
text(0.0005, 0.85, paste( "Occurrence: ", round(occ_frac, digits=4) ))
dev.off()

png( file = paste( dir, dataname, '_SFCurve_MOPS.png', sep = "" ) )
plot(fdr_mops[list_mops], recall_mops[list_mops], log="x",
     main=paste(dataname, "Sensitivity vs. FDR for MOPS model", sep="\n\n"), 
     xlab="FDR", ylab="Sensitivity", xlim=c(0.0001,0.5), ylim=c(0,1),
     type='l', lwd=2.5,
     col="deepskyblue")
if(right_mops == left_mops){
  ausfc_mops = 0 
} else {
  ausfc_mops = sum(diff(log10(fdr_mops[list_mops]))*rollmean(recall_mops[list_mops],2)) / sum_area
}

text(0.0005, 0.95, paste( "AUSFC: ", round(ausfc_mops, digits=4) ))
text(0.0005, 0.85, paste( "Mult. Occ.: ", round(occ_mult, digits=4) ))
dev.off()

# write paramaters into a file
ofile = paste( dir, dataname, '.rankscore', sep = "" )
writeLines(c(
             "fraction occurrence: ", occ_frac, 
             
             "multiple occurrence: ", occ_mult,
            
             "eta0 (MOPS): ", eta0_mops,
             
             "eta0 (ZOOPS): ", eta0_zoops,
             
             "AUSFC (MOPS): ", ausfc_mops,
             
             "AUSFC (ZOOPS): ", ausfc_zoops,
            
             "Ranking score (MOPS): ", ausfc_mops*occ_mult,
             
             "Ranking score (ZOOPS): ", ausfc_zoops*occ_frac
             ), ofile)

