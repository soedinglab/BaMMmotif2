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
# ./plotAUSFC_benchmark_fdrtool.R PATH_TO_zoops.stats_FILE BASENAME_OF_THE_FILE PATH_OF_OUTPUT_FILE
# ./plotAUSFC_benchmark_fdrtool.R /home/bamm_result/ JunD /home/bamm_result/output.txt

#-----------------------------
#
# R sources
#
#-----------------------------

# load "fdrtool" library for calculating FDR and recall
library( fdrtool )

# load "zoo" package for calculating AUPRC
library( zoo )

# load "argparse" library for parsing arguments
library( argparse )

#-----------------------------
#
# Parameter setting
#
#-----------------------------
parser <- ArgumentParser(description='benchmark motif finder')
parser$add_argument('target_directory', help='directory that contains the target file')
parser$add_argument('prefix', help='prefix of the target file')
parser$add_argument('output_file', help='file to print the benchmark results')

args <- parser$parse_args()

# default args for ArgumentParser()$parse_args are commandArgs(TRUE)
# which is what you'd want for an Rscript but not for interactive use

# get the directory of the p-value files
dir <- args$target_directory
prefix <- args$prefix
output <- args$output_file

results = c(paste(c("prefix", "motifNumber", "ausfc_mops*occ_mult", "ausfc_zoops*occ_frac", "occ_mult", "occ_frac"), collapse="\t"))
for (f in Sys.glob(paste(c(dir, "/", prefix, "_motif_*", ".zoops.pvalues"), collapse=""))) {
  motifNumber <- sub(paste(c(dir, "/", prefix, "_motif_"), collapse=""), "", f)
  motifNumber <- sub(".zoops.pvalues", "", motifNumber)

  ########### for zoops model ###########
  # read in p-values from the files
  pvalues_zoops <- unlist(read.table(paste(c(dir, "/", prefix, "_motif_", motifNumber, ".zoops.pvalues"), collapse="")))
  
  # avoid the rounding errors when p-value = 0 or p-value > 1
  for(i in seq(1, length(pvalues_zoops))){
    if( pvalues_zoops[i] > 1 ){
      pvalues_zoops[i] = 1
    }
  }
  
  # estimate False Discovery Rates for Diverse Test Statistics using fdrtool
  result_zoops = fdrtool( pvalues_zoops, statistic="pvalue",
                          plot=FALSE, color.figure=TRUE, verbose=TRUE )
 
  # get the global fdr values and estimate of the weight eta0 of
  # the null component
  fdr_zoops <- result_zoops$qval
  eta0_zoops <- result_zoops$param[3]
  
  # avoid eta0 = 1
  if(eta0_zoops==1){
    eta0_zoops=0.9999999
  }
  
  # calculate recall
  len_zoops = length(fdr_zoops)
  list_zoops <- seq(1, len_zoops )
  recall_zoops <- ( 1 - fdr_zoops ) * list_zoops / ( 1 - eta0_zoops ) / len_zoops
  
  # modify the SF curve:
  # reset recall to 1 when it is larger than 1
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
  list_zoops <- seq(1, len_zoops+1)
  
  # limit the frame of the curve
  left_zoops = 1
  right_zoops = len_zoops+1
  
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
  
  sum_area = log10(0.5) + 4
  
  # calculate the occurrence fraction of motifs
  occ_frac = 1 - eta0_zoops
  
  # compute  the area under the sensitivity-FDR curve(AUSFC):
  if(right_zoops == left_zoops){
    ausfc_zoops = 0
  } else {
    ausfc_zoops = sum(diff(log10(fdr_zoops[list_zoops]))*rollmean(recall_zoops[list_zoops],2)) / sum_area
  }
  # plot fdr vs. recall curve
  png( file = paste( dir, "/", prefix, "_motif_", motifNumber, '_SFCurve_ZOOPS.png', sep = "" ) )
  plot(fdr_zoops[list_zoops], recall_zoops[list_zoops], log="x",
       main=paste(prefix, "_motif_", motifNumber, " Sensitivity vs. FDR for ZOOPS model", sep=""),
       xlab="FDR", ylab="Sensitivity", xlim=c(0.0001,0.5), ylim=c(0,1),
       type='l', lwd=2.5,
       col="deepskyblue")
  text(0.0005, 0.95, paste( "AUSFC: ", round(ausfc_zoops, digits=4) ))
  text(0.0005, 0.85, paste( "Frac. Occ: ", round(occ_frac, digits=4) ))
  dev.off()
  
  
  ########### for mops model ###########
  # read in p-values from the files
  pvalues_mops <- unlist(read.table(paste(c(dir, "/", prefix, "_motif_", motifNumber, ".mops.pvalues"), collapse="") ))

  # avoid the rounding errors when p-value = 0 or p-value > 1
  for(i in seq(1, length(pvalues_mops))){
    if( pvalues_mops[i] > 1 ){
      pvalues_mops[i] = 1
    }
  }
  
  # estimate False Discovery Rates for Diverse Test Statistics using fdrtool
  result_mops = fdrtool( pvalues_mops, statistic="pvalue",
                         plot=FALSE, color.figure=TRUE, verbose=TRUE )

  # get the global fdr values and estimate of the weight eta0 of
  # the null component
  fdr_mops <- result_mops$qval
  eta0_mops <- result_mops$param[3]

  # avoid eta0 = 1
  if(eta0_mops==1){
    eta0_mops=0.9999999
  }

  # calculate recall
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

  # extend the SF curve till FDR reaches 0.5
  recall_mops <- append(recall_mops, 1)
  fdr_mops <- append(fdr_mops, 0.5)
  list_mops <- seq(1, len_mops+1)

  # limit the frame of the curve
  left_mops = 1
  right_mops = len_mops+1
  
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

  # calculate motif occurrences per sequence
  occ_mult = (1 - eta0_mops) * len_mops / len_zoops

  # compute  the area under the sensitivity-FDR curve(AUSFC):
  if(right_mops == left_mops){
    ausfc_mops = 0
  } else {
    ausfc_mops = sum(diff(log10(fdr_mops[list_mops]))*rollmean(recall_mops[list_mops],2)) / sum_area
  }
  
  # plot fdr vs. recall curve
  png( file = paste( dir, "/", prefix, "_motif_", motifNumber, '_SFCurve_MOPS.png', sep = "" ) )
  plot(fdr_mops[list_mops], recall_mops[list_mops], log="x",
       main=paste(prefix, "_motif_", motifNumber, " Sensitivity vs. FDR for MOPS model", sep=""),
       xlab="FDR", ylab="Sensitivity", xlim=c(0.0001,0.5), ylim=c(0,1),
       type='l', lwd=2.5,
       col="deepskyblue")
  text(0.0005, 0.95, paste( "AUSFC: ", round(ausfc_mops, digits=4) ))
  text(0.0005, 0.85, paste( "Mult. Occ.: ", round(occ_mult, digits=4) ))
  dev.off()

  resultString = paste(c(prefix, motifNumber, ausfc_mops*occ_mult, ausfc_zoops*occ_frac, occ_mult, occ_frac), collapse="\t")
  print(resultString)
  results = c(results, resultString)
}

outConn <- file(output)
writeLines(results, outConn)
close(outConn)
