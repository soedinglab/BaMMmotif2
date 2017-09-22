#!/usr/bin/env Rscript
#-----------------------------------------------------------------
#
# Documentation
#
#-----------------------------------------------------------------
# based on the output of BaMM FDR function:
# 1. calculate false discovery rate(FDR), and sensitivity
#    using the p-values for log odds ratios of positive sequences
#    and Grenander fitting in the "fdrtool" library;
# 2. calculate area-under-FDR-sensitivity-curve (AUSFC);
# 3. plot FDR-sensitivity curve using fitted FDR and sensitivity scores;
# 4. plot partial ROC using raw TPs and FPs;
# 5. plot precision-recall curve using raw FDR and recall scores;
# 6. ausfc score, pAUC score and auprc score are saved in .bmscore file.

# command for executing this script:
# ./evaluateBaMM.R INPUT_DIR FILE_PREFIX [options]

# ./evaluateBaMM.R /home/bamm_result/ JunD_motif_1 [--SFC 1 ...]

#-----------------------------
#
# R sources
#
#-----------------------------

# load "zoo" library for calculating area-under-the-curve
library( zoo )
# load "argparse" library for parsing arguments
library( argparse )
# load "fdrtool" library for calculating FDR 
library( fdrtool )
# load "LSD" library for plotting the curves
library( LSD )

#-----------------------------
#
# Parameter setting
#
#-----------------------------
parser <- ArgumentParser(description="evaluate motifs optimized by BaMM")
# positional arguments
parser$add_argument('target_directory', help="directory that contains the target file")
parser$add_argument('prefix', help="prefix of the target file")
# optional arguments
parser$add_argument("--fdrtool", type="logical", default=FALSE, help="flag for print out analysis fromfdrtool" )
parser$add_argument("--SFC", type="logical", default=FALSE, help="flag for print out sensitivity-fdr curve" )
parser$add_argument("--ROC5", type="logical", default=FALSE, help="flag for print out partial ROC curve" )
parser$add_argument("--PRC", type="logical", default=FALSE, help="flag for print out precision-recall curve" )

args <- parser$parse_args()

# get the directory of the p-value files
dir 	<- args$target_directory
prefix 	<- args$prefix
# flags for printing the curve plots
print_FDRtool	= args$fdrtool
print_SFcurve 	= args$SFC
print_ROC5      = args$ROC5
print_PRcurve 	= args$PRC

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


#-----------------
#
# main functions
#
#-----------------
results = c()
#resultTitle = paste0(c("motif_name", "ausfc_score", "auc5_score", "auprc_score"), collapse="\t")
#results = c(results, resultTitle)

for (f in Sys.glob(paste(c(dir, "/", prefix, "*", ".zoops.stats"), collapse=""))) {

    # read in p-values from file
    first_row   <- read.table(f, nrows =1 )
    stats       <- read.table(f, skip=1 )
    pvalues     <- stats$V5
    mfold       <- as.numeric(first_row[6])

    # avoid the rounding errors when p-value = 0 or p-value > 1
    for(i in seq(1, length(pvalues))){
        if( pvalues[i] > 1 ){
            pvalues[i] = 1
        }
    }

    # estimate False Discovery Rates for Diverse Test Statistics
    eta0set = mfold / ( 1+mfold )
    if( print_FDRtool ){
        pdf( file = paste(dir, prefix, '_FDRstat.pdf', sep = "" ) )
        result = fdrtool( pvalues, statistic="pvalue",
                          plot=TRUE, color.figure=TRUE, verbose=TRUE, eta0set=eta0set )
        dev.off()
    } else {
        result = fdrtool( pvalues, statistic="pvalue", plot=FALSE, eta0set=eta0set )
    }

    # get the global fdr values and estimate of the weight eta0 of
    # the null component
    fdr_m	<- result$qval
    eta0 	<- result$param[3]

    #calculate FDR
    fdr 	<- 1 / ( 1 + mfold * ( 1 / fdr_m - 1 ) )

    # calculate recall
    len 	= length(fdr)
    list 	<- seq(1, len)
    recall  <- ( 1 - fdr ) * list / ( 1 - eta0 ) / len

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
    recall  <- append(recall, 1)
    fdr     <- append(fdr, 0.5)
    range   <- seq(1, len+1)

    # limit the frame of the curve to FDR(0.01-0.5)
    left    = 1         # index for the FDR value on the left border
    right   = len+1     # index for the FDR value on the right border
    l_range = 0.01	    # left range for FDR
    r_range = 0.5		# right range for FDR

    for(i in range){
        if( fdr[i] >= l_range ){
            left = i
            recall[i] = 0
            break
        }
    }

    for(i in range){
        if( fdr[i] > r_range ){
            right = i
            break
        }
    }

    range <- seq(left, right)

    sum_area = log10(r_range)-log10(l_range)

    # compute the area under the sensitivity-FDR curve (AUSFC):
    if(right == left){
        ausfc = 0
    } else {
        ausfc = sum(diff(log10(fdr[range]))*rollmean(recall[range],2)) / sum_area
    }

    #####################################################################
    #
    # plot the Sensitivity-FDR curve, with FDR in the range of (0.01~0.5)
    #
    ######################################################################
    if( print_SFcurve ){
        # plot fdr vs. recall curve. Note that x-axis is in log scale
        picname <- paste0( dir, '/', prefix, '_SFC.jpeg' )
        jpeg( filename = picname, width = 800, height = 800, quality = 100 )
        par(oma=c(0,0,0,0), mar=c(6,6.5,5,1))
        plot(fdr[range], recall[range],
            log="x",
            main="Motif Performance",
            xlim=c(l_range,r_range), ylim=c(0,1),
            xlab="", ylab="",
            type="l", lwd=7.5,
            col="orange",
            axes = FALSE, cex.main=3.0
        )
        mtext("False discovery rate", side=1, line=4.5, cex = 3.5)
        mtext("Sensitivity", side=2, line=4, cex = 3.5)
        axis(1, at=c(0.01,0.05,0.5),labels = c(0.01,0.05,0.5),tick = FALSE, cex.axis=3.0, line=1)
        axis(2, at=c(0,0.5,1),labels = c(0,0.5,1),tick = FALSE, cex.axis=3.0)
        polygon(c(l_range,fdr[range],r_range), c(0, recall[range],0), col = convertcolor("orange",30), border = NA)
        text(x = 0.1,y = 0.1,labels = paste0("AUSFC = ", round(ausfc,digits=3)), cex = 3.5)
        box(lwd=2.5)
        invisible(dev.off())
    }

    ########################################################
    #
    #   calculate AUC5 under the ROC curve with FDR~(0,0.05)
    #
    ########################################################
    TP  <- stats$V1
    FP  <- stats$V2
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
        picname <- paste0( dir, '/', prefix, '_pROC.jpeg' )
        jpeg( filename = picname, width = 800, height = 800, quality = 100 )
        par(oma=c(0,0,0,0), mar=c(6,6.5,5,1))
        plot(FPR[1:rbound], TPR[1:rbound],
            main="Partial ROC",
            xlim=range(0,0.05), ylim=range(0,1),
            xlab="", ylab="",
            type="l", lwd=7.5,
            col="darkgreen",
            axes = FALSE, cex.main=3.0
        )
        mtext("False positive rate", side=1, line=4.5, cex = 3.5)
        mtext("True positive rate", side=2, line=4, cex = 3.5)
        axis(1, at=c(0,0.025,0.05),labels = c(0,0.025,0.05),tick = FALSE, cex.axis=3.0, line=1)
        axis(2, at=c(0,0.5,1),labels = c(0,0.5,1),tick = FALSE, cex.axis=3.0)
        polygon(c(0,FPR[1:rbound],0.05), c(0,TPR[1:rbound],0), col = convertcolor("darkgreen",30), border = NA)
        text(x = 0.025,y = 0.1,labels = paste0("pAUC = ", round(auc5,digits=3)), cex = 3.5)
        box(lwd=2.5)
        invisible(dev.off())
    }

    ###################################################################
    #
    # plot the Precision-Recall curve using the raw values from BaMM
    #
    ###################################################################
    raw_fdr         <- stats$V3
    raw_recall      <- stats$V4
    # solve numeric issue
    # reset recall to 1 when it is larger than 1
    for(i in seq(1,length(raw_recall))){
        if( raw_recall[i] > 1 ){
            for(rest in i:length(raw_recall)){
                raw_recall[rest] = 1
            }
            break
        }
    }

    raw_precision   = 1 - raw_fdr
    auprc           = sum(diff(raw_recall)*rollmean(raw_precision,2))
    # plot the raw precision-recall curve
    if( print_PRcurve ){
        picname <- paste0( dir, '/', prefix, '_PRC.jpeg' )
        jpeg( filename = picname, width = 800, height = 800, quality = 100 )
        par(oma=c(0,0,0,0), mar=c(6,6.5,5,1))
        plot(raw_recall, raw_precision,
            main="Motif Performance",
            xlim=range(0,1), ylim=range(0,1),
            xlab="", ylab="",
            type="l", lwd=7.5,
            col="darkblue",
            axes = FALSE, cex.main=3.0
        )
        mtext("Recall", side=1, line=4.5, cex = 3.5)
        mtext("Precision", side=2, line=4, cex = 3.5)
        axis(1, at=c(0,0.5,1),labels = c(0,0.5,1),tick = FALSE, cex.axis=3.0, line=1)
        axis(2, at=c(0,0.5,1),labels = c(0,0.5,1),tick = FALSE, cex.axis=3.0)
        polygon(c(0,raw_recall,1), c(0,raw_precision,0), col = convertcolor("darkblue",30), border = NA)
        text(x = 0.5,y = 0.3,labels = paste0("AUC = ", round(auprc,digits=3)), cex = 3.5)
        box(lwd=2.5)
        invisible(dev.off())
    }

    resultString = paste0(c(prefix, ausfc, auc5, auprc), collapse="\t")
    #print( resultString )
    results = c(results, resultString)
}
outConn <- file(paste(dir, '/', prefix, ".bmscore", sep = "" ))
writeLines(results, outConn)
close(outConn)
