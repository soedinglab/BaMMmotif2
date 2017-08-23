#!/usr/bin/env Rscript
#-----------------------------------------------------------------
#
# Documentation
#
#-----------------------------------------------------------------

# calculate FDR, precision, recall and occurrences for both MOPS
# and ZOOPS models by taking p-values of log odds scores froom
# positive sequences and using "fdrtool" library
# plot FDR vs. recall(sensitivity) curve
# calculate the area under the sensitivity-FDR curve (AUSFC)
# and plot the true-positive-rate(TPR) vs. false-positive-rate(FPR)
# curve and calculate the partial AUC from this curve;
# results are saved in a .bmscore file.

# examples for running this script:
# ./plotAUSFC_benchmark_fdrtool.R PATH_TO_zoops.stats_FILE BASENAME_OF_THE_FILE PATH_TO_OUTPUT_FILE
# ./plotAUSFC_benchmark_fdrtool.R /home/bamm_result/ JunD_motif_1 /home/bamm_result/ausfc.txt

#-----------------------------
#
# R sources
#
#-----------------------------

# load "zoo" library for calculating AUPRC
library( zoo )
# load "argparse" library for parsing arguments
library( argparse )

###########################
## For this script we need a slightly modified version of the fdrtool function.
## Unfortunately, R does not seem to provide an easy way to source an R file
## in the same directory. We sincerely apologize for copying the source here.
## Please forgive us. We won't do it again. We promise.

### Modified version by Wanwan Ge (wanwan.ge@mpibpc.mpg.de)
###
### Copyright 2017 Wanwan Ge
###
### fdrtool.R  (2013-09-15)
###
###    Estimate (Local) False Discovery Rates For Diverse Test Statistics
###
###
### Copyright 2007-13 Korbinian Strimmer
###
### This file is part of the `fdrtool' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

library("fdrtool")

fdrtool = function(x,
  statistic=c("normal", "correlation", "pvalue"),
  #statistic=c("normal", "correlation", "pvalue", "studentt"),
  plot=TRUE, color.figure=TRUE, verbose=FALSE,
  cutoff.method=c("fndr", "pct0", "locfdr"),
  pct0=0.75,
  eta0set=0.9091)
{
  statistic = match.arg(statistic)
  cutoff.method = match.arg(cutoff.method)

  if ( is.vector(x) == FALSE )
  	stop("input test statistics must be given as a vector!")

  if ( length(x) < 200 ) warning("There may be too few input test statistics for reliable FDR calculations!")
  if (statistic=="pvalue")
  {
    if (max(x) > 1 | min(x) < 0)
      stop("input p-values must all be in the range 0 to 1!")
  }

#### step 1 ####

  if(verbose) cat("Step 1... determine cutoff point\n")

  # determine cutoff point for censoring

  if (cutoff.method=="pct0")
  {
    # use specified quantile

    if(statistic=="pvalue") x0 = quantile(x, probs=1-pct0)
    else x0 = quantile(abs(x), probs=pct0)
  }
  else if ( cutoff.method=="locfdr" & (statistic=="normal" | statistic=="correlation") )
  {
    # use same procedure as in locfdr R package (due to Brit Katzen-Turnbull)

    if(statistic=="normal") z = x
    if(statistic=="correlation") z = atanh(x)

    iqr = as.double(diff(quantile(z, probs=c(.25, .75))))
    sdhat = iqr/(2*qnorm(.75))
    N = length(z)
    # b = 3.55-.44*log(N, 10)                               # locfdr 1.1-3
    b = ifelse(N > 500000, 1,  4.3 * exp(-0.26*log(N,10)) ) # locfdr 1.1-6
    z0 = b*sdhat

    if(statistic=="normal") x0 = z0
    if(statistic=="correlation") x0 = tanh(z0)
  }
  else
  {
    if(cutoff.method=="locfdr")
    warning("cutoff.method=\"locfdr\" only available for normal and correlation statistic.")

    # control false nondiscovery rate

    x0 = fndr.cutoff(x, statistic)
  }

#### step 2 ####

  if(verbose) cat("Step 2... estimate parameters of null distribution and eta0\n")

  cf.out <- censored.fit(x=x, cutoff=x0, statistic=statistic)
# cf.out looks as follows for p-values
#     cutoff N0      eta0    eta0.var
#[1,]   0.96 64 0.3730473 0.002141996
# the other test statistics have two more colums containing the scale parameter

  if (statistic=="pvalue")
    scale.param = NULL
  else
    scale.param <- cf.out[1,5] # variance parameter

  cf.out[1,3] <- eta0set		# <---- Here is what is changed by the modified version
  eta0 = cf.out[1,3]


#### step 2 ####

  if(verbose) cat("Step 3... compute p-values and estimate empirical PDF/CDF\n")

  nm = fdrtool:::get.nullmodel(statistic)
  pval = nm$get.pval(x, scale.param)

  # determine cumulative empirical distribution function (pvalues)
  ee <- fdrtool:::ecdf.pval(pval, eta0=eta0)

  g.pval <- grenander(ee)

  #cat("DEBUG: Grenander eta0=", g.pval$f.knots[length(g.pval$f.knots)], "\n")
  #cat("DEBUG: estimated eta0=", eta0 , "\n\n")

  # mixture density and CDF
  f.pval = approxfun( g.pval$x.knots,  g.pval$f.knots, method="constant", rule=2)
  f0.pval = function(x) return( ifelse(x > 1 | x < 0, 0, rep(1, length(x))) )

  F.pval = approxfun( g.pval$x.knots,  g.pval$F.knots, method="linear",
           yleft=0, yright=g.pval$F.knots[length(g.pval$F.knots)])
  F0.pval = function(x) return( ifelse(x > 1, 1, ifelse(x < 0, 0, x )) )

  #fdr.pval = function(p) pmin( eta0   / f.pval(p), 1) # eta0*f0/ f
  fdr.pval = function(p)
  {
    p[ p == .Machine$double.eps ] = 0
    pmin( eta0   / f.pval(p), 1) # eta0*f0/ f
  }

  Fdr.pval = function(p) pmin( eta0*p / F.pval(p), 1) # eta0*F0/ F


#### step 4 ####

  if(verbose) cat("Step 4... compute q-values and local fdr\n")

  qval <- Fdr.pval(pval)
  lfdr <- fdr.pval(pval)

#### return results ####

  result = list(pval=pval, qval=qval, lfdr=lfdr,
             statistic=statistic, param=cf.out)

  if (plot)
  {
    if(verbose) cat("Step 5... prepare for plotting\n")

    ##############
    # zeta > 0 in the following

    if(statistic=="pvalue")
    {
      f0 <- function(zeta) return( nm$f0(zeta, scale.param) )
      F0 <- function(zeta) return( nm$F0(zeta, scale.param) )
      get.pval <- function(zeta) return( nm$get.pval(1-zeta, scale.param) )
      x0 = 1-x0
    }
    else
    {
      f0 <- function(zeta) return( 2*nm$f0(zeta, scale.param)  )
      F0 <- function(zeta) return( 2*nm$F0(zeta, scale.param)-1  )
      get.pval <- function(zeta) return( nm$get.pval(zeta, scale.param) )
    }

    fdr = function(zeta)  fdr.pval(get.pval(zeta))
    Fdr = function(zeta)  Fdr.pval(get.pval(zeta))

    F   = function(zeta)  1-eta0*get.pval(zeta)/Fdr(zeta)
    FA  = function(zeta)  (F(zeta)-eta0*F0(zeta))/(1-eta0)

    f   = function(zeta)  eta0*(f0(zeta))/fdr(zeta)
    fA  = function(zeta)  (f(zeta)-eta0*f0(zeta))/(1-eta0)

    ##############

    ax = abs(x)
    if (statistic=="pvalue") ax = 1-ax  # reverse p-val plot

    xxx = seq(0, max(ax), length.out=500)

    ll = pvt.plotlabels(statistic, scale.param, eta0)

    par(mfrow=c(3,1))

    if (color.figure)
      cols = c(2,4) # colors for f0,F0 and fA,FA
    else
      cols = c(1,1)

    hist(ax, freq=FALSE, bre=50,
      main=ll$main, xlab=ll$xlab, cex.main=1.8)
    lines(xxx, eta0*f0(xxx), col=cols[1], lwd=2, lty=3 )
    lines(xxx, (1-eta0)*fA(xxx), col=cols[2], lwd=2 )
    if (statistic=="pvalue")
      pos1 = "topleft" else pos1="topright"
    legend(pos1,
      c("Mixture", "Null Component", "Alternative Component"),
      lwd=c(1, 2, 2), col=c(1,cols), lty=c(1,3,1), bty="n", cex=1.5)

    plot(xxx, F(xxx), lwd=1, type="l", ylim=c(0,1),
      main="Density (first row) and Distribution Function (second row)",
      xlab=ll$xlab, ylab="CDF", cex.main=1.5)
    lines(xxx, eta0*F0(xxx), col=cols[1], lwd=2, lty=3)
    lines(xxx, (1-eta0)*FA(xxx), col=cols[2], lwd=2)

    # DEBUG show cutoff in green line
    #lines(c(x0,x0),c(0,1), col=3)

    plot(xxx, Fdr(xxx), type="l", lwd=2, ylim=c(0,1),
      main="(Local) False Discovery Rate", ylab="Fdr and fdr",
      xlab=ll$xlab, lty=3, cex.main=1.5)
    lines(xxx, fdr(xxx), lwd=2)
    if (eta0 > 0.98)
      pos2 = "bottomleft" else pos2="topright"
    legend(pos2,
      c("fdr (density-based)", "Fdr (tail area-based)"),
      lwd=c(2,2), lty=c(1,3), bty="n", cex=1.5)

    # DEBUG show cutoff in green line
    #lines(c(x0,x0),c(0,1), col=3)

    par(mfrow=c(1,1))

    rm(ax)
  }

  if(verbose) cat("\n")

  return(result)
}

#####

## create labels for plots
pvt.plotlabels <- function(statistic, scale.param, eta0)
{
   if (statistic=="pvalue")
   {
     main = paste("Type of Statistic: p-Value (eta0 = ", round(eta0, 4), ")", sep="")
     xlab ="1-pval"
   }

   if (statistic=="studentt")
   {
     df = scale.param
     main = paste("Type of Statistic: t-Score (df = ", round(df,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(t)"
   }

   if (statistic=="normal")
   {
     sd = scale.param
     main = paste("Type of Statistic: z-Score  (sd = ", round(sd,3),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(z)"
   }

   if (statistic=="correlation")
   {
     kappa =scale.param
     main = paste("Type of Statistic: Correlation (kappa = ", round(kappa,1),
                       ", eta0 = ", round(eta0, 4), ")", sep="")
     xlab = "abs(r)"
   }

   return(list(main=main, xlab=xlab))
}
## EOF fdrtool.R

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

# flags for printing the curve plots
print_PRcurve = FALSE
print_FDRcurve = FALSE
print_SFcurve = FALSE
print_ROC5 = FALSE

results = c()
for (f in Sys.glob(paste(c(dir, "/", prefix, "*", ".zoops.stats"), collapse=""))) {
  motifNumber <- sub(paste(c(dir, "/", prefix, "_motif_"), collapse=""), "", f)
  motifNumber <- sub(".zoops.stats", "", motifNumber)

  # read in p-values from file
  first_row <- read.table(f, nrows = 1 )
  stats <- read.table(f, skip=1 )
  pvalues <- stats$V5
  mfold <- as.numeric( first_row[6] )

  # avoid the rounding errors when p-value = 0 or p-value > 1
  for(i in seq(1, length(pvalues))){
    if( pvalues[i] > 1 ){
      pvalues[i] = 1
    }
  }

  # estimate False Discovery Rates for Diverse Test Statistics
  eta0set = mfold / (1+mfold)
  if( print_FDRcurve ){
    pdf( file = paste(dir, '/', prefix, "_motif_", motifNumber, '_FDRstat.pdf', sep = "" ) )
    result = fdrtool( pvalues, statistic="pvalue",
                      plot=TRUE, color.figure=TRUE, verbose=TRUE, eta0set=eta0set )
    dev.off()
  } else {
    result = fdrtool( pvalues, statistic="pvalue", plot=FALSE, eta0set=eta0set )
  }

  # get the global fdr values and estimate of the weight eta0 of
  # the null component
  fdr_m <- result$qval
  eta0 <- result$param[3]
  fdr <- 1 / ( 1 + mfold * ( 1 / fdr_m - 1 ) )

  # calculate recall
  len = length(fdr)
  list <- seq(1, len )
  recall <- ( 1 - fdr ) * list / ( 1 - eta0 ) / len

  # calculate precision 
  raw_fdr <- stats$V3
  raw_recall <- stats$V4
  precision = 1 - raw_fdr
  if( print_PRcurve ){
    pdf( file = paste(dir, '/', prefix, "_motif_", motifNumber, '_PRcurve.pdf', sep = "" ) )
    plot(raw_recall, precision,
         main=paste(prefix, "_motif_", motifNumber, "\n Precision-Recall curve", sep = ""),
         xlab="Recall", ylab="Precision", xlim=c(0,1), ylim=c(0,1),
         type='l', lwd=2.5, col="green")
    # compute the area under the partical TP-FP curve(AUC5):
    auprc = sum(diff(precision)*rollmean(raw_recall,2))
    
    # write the AUC on the plot with the precision of 4 digits:
    text(0.1, 0.9, paste( "AUPRC: ", round(auprc, digits=4) ))
    dev.off()
  }
  
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
    if( fdr[i] >= 0.01 ){
      left = i
      recall[i] = 0
      break
    }
  }

  for(i in range){
    if(fdr[i] > 0.5 ){
      right = i
      break
    }
  }

  range <- seq(left, right)

  sum_area = log10(0.5) + 2
  # compute  the area under the sensitivity-FDR curve(AUSFC):
  if(right == left){
    ausfc = 0
  } else {
    ausfc = sum(diff(log10(fdr[range]))*rollmean(recall[range],2)) / sum_area
  }
  
  if( print_SFcurve ){
    # plot fdr vs. recall curve
    pdf( file = paste(dir, '/', prefix, "_motif_", motifNumber, '_SFCurve.pdf', sep = "" ) )
    plot(fdr[range], recall[range], log="x",
         main=paste(prefix, "_motif_", motifNumber, "\nSensitivity vs. FDR", sep=""),
         xlab="FDR", ylab="Sensitivity", xlim=c(0.01,0.5), ylim=c(0,1),
         type='l', lwd=2.5,
         col="deepskyblue")
    # fill the area under the FDR-recall curve 
    #  polygon(c(min(fdr), fdr, max(fdr)), c(min(recall),recall,1), col="gray", border="gray")
    #  polygon(c(min(fdr), range, 0.5), c(min(recall),recall,1), col="gray", border="gray")
    # write the AUSFC on the plot with the precision of 4 digits:
    text(0.35, 0.05, paste( "AUSFC: ", round(ausfc, digits=4) ))
    dev.off()
  }

  # calculate AUC5 under the ROC curve
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
  # compute the area under the partical TP-FP curve(AUC5):
  auc5 = sum(diff(FPR[1:rbound])*rollmean(TPR[1:rbound],2)) / 0.05

  if( print_ROC5 ){
    pdf( file = paste(dir, '/', prefix, "_motif_", motifNumber, '_pROC.pdf', sep = "" ) )
  
    plot(FPR[1:rbound], TPR[1:rbound],
           main=paste(prefix, "_motif_", motifNumber, " FPR vs. TPR ", sep=""),
           xlab="FPR", ylab="TPR", xlim=c(0,0.05), ylim=c(0,1),
           type='l', lwd=2.5, col="deepskyblue")
    # write the AUC on the plot with the precision of 4 digits:
    text(0.03, 0.1, paste( "pAUC: ", round(auc5, digits=4) ))
    dev.off()
  }
  resultString = paste(c(prefix, motifNumber, ausfc, auc5), collapse="\t")
  #print(resultString)
  results = c(results, resultString)
}

outConn <- file(output)
writeLines(results, outConn)
close(outConn)
