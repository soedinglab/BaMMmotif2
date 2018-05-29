#!/usr/bin/env Rscript
#-----------------------------------------------------------------
#
# Documentation
#
#-----------------------------------------------------------------
# based on the output of BaMM FDR function (.stats file):
# 1. p-value statistics from FDR for both motif itself and motif tested on the dataset
# 2. plot Recall-TP/FP ratio curves for both motif itself and motif tested on the dataset
# 3. output the AvRec score and motif occurrence in a .bmscore file

# command for executing this script:
# ./plotPvalStat.R INPUT_DIR FILE_PREFIX [options]

# ./plotPvalStat.R /home/bamm_result/ JunD [--web 1] [--plots 1]

#-----------------------------
#
# R sources
#
#-----------------------------

# load "zoo" library for calculating area-under-the-curve
suppressMessages(library( zoo ))
# load "argparse" library for parsing arguments
suppressMessages(library( argparse ))
# load "fdrtool" library for calculating FDR
suppressMessages(library( fdrtool ))
# load "LSD" library for plotting the curves
suppressMessages(library( LSD ))

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
parser$add_argument("--plots", type="logical", default=FALSE, help="flag for printing out plots" )
parser$add_argument("--web", type="logical", default=FALSE, help="flag for printing out AvRec score on the screen" )

# parse the arguments
args    <- parser$parse_args()

# interpret the arguments
dir 	<- args$target_directory
prefix 	<- args$prefix

# flag for verbose output for web usage
web     <- args$web

# flag for plots
plots   <- args$plots

# preset variables
ofile   <- paste(dir, '/', prefix, ".bmscore", sep = "" )

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
eta0set=NULL)
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

    if( !is.null(eta0set) )  cf.out[1,3] <- eta0set		# <---- Here is what is changed by the modified version
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


#-----------------------------------------
#
# plot functions used for motif evaluation
#
#-----------------------------------------

# plot p-value statistics using fdrtool
plotPvalStat = function(pvalues, filename, eta0, data_eta0, rerank){

    picname = paste0(filename,".png")
    png( filename = picname, width = 800, height = 800 )
    histRes <- hist(pvalues, plot=FALSE, breaks=30)
    xvals <- histRes$breaks
    yvals <- histRes$density

    xvals <- c(xvals, 0)
    yvals <- c(0, yvals, 0)

    if(rerank){
        mainname = "Motif P-value statistics"
    } else{
        mainname = "Dataset P-value statistics"
    }

    par(oma=c(0,0,0,0), mar=c(6,6.5,5,1))

    rbound = max(xvals)

    plot(xvals,yvals,
    type="S",
    main=mainname,
    xlab="", ylab="",
    lwd=3, axes=FALSE, frame.plot=TRUE, cex.main = 3.0, cex.axis=2.5)

    mtext("P-values", side=1, line=4.5, cex = 3)
    mtext("Density", side=2, line=4, cex = 3)

    #axis(1,tick =FALSE, cex.axis=2.5, line=1)
    axis(1, at=c(0,0.5,1), labels = c(0,0.5,1), tick =FALSE, cex.axis=2.5, line=1)

    # mark the negative regions from background sequences
    rect(0, 0, rbound, data_eta0,
        col="grey", angle=45, density=4, lty=NULL, lwd=2, xpd=FALSE)
        lines(xvals,yvals,type="S")

    cutoff=0.1      # cutoff for p-values
    abline(v=cutoff, col="red", lwd=4)

    if(rerank){
        abline(h=eta0, col="orange", lwd=5, lty=2)
    }

    # color FP region
    rect(0, 0, cutoff, eta0,
        col=rgb(1, 0, 0, 0.1),
        border=par("fg"), lty=NULL, lwd=0, xpd=FALSE)

    mask = xvals <= cutoff

    xpoly = xvals[mask]
    ypoly = yvals[mask]

    xpoly = rep(xpoly[2:length(xpoly)], each=2)
    ypoly = c(rep(ypoly[2:length(ypoly)], each=2), rep(ypoly[length(ypoly)], each=2))
    xpoly = c(0, xpoly[1:length(xpoly) - 1], cutoff, cutoff)
    ypoly = pmax(ypoly, eta0)

    # color TP region
    polygon(
    c(xpoly, rev(xpoly)), c(ypoly, rev(rep(eta0, length(ypoly)))),
    col =rgb(0, 1, 0, 0.1),
    border = NA
    )

    text_cex = 2
    font = 2
    v_spacer = 0.05
    #h_spacer = max(ypoly) * 0.04  # scale the spacer due to the range of y-axis
    h_spacer = text_cex * 0.1

    text(cutoff - v_spacer, eta0 + h_spacer, "TP", col="darkgreen", font=font, cex=text_cex)
    text(cutoff + v_spacer, max(0, eta0 - h_spacer), "TN", col="black", font=font, cex=text_cex)
    text(cutoff - v_spacer, max(0, eta0 - h_spacer), "FP", col="darkred", font=font, cex=text_cex)
    text(cutoff + v_spacer, eta0 + h_spacer, "FN", col="black", font=font, cex=text_cex)

    text(rbound/2, data_eta0/2, "background sequences", font=font, cex=text_cex, col="gray30")

    box(lwd=2.5)

    invisible(dev.off())
}


# plot TP/FP ratio vs. recall curve
plotRRC = function(picname, recall, TFR, rerank){

    # set the upper and lower region for the y-axis
    y_upper = 100
    y_lower = 1

    # compute the area under the RRC curve (AvRec):
    sum_area = log10(y_upper)-log10(y_lower)    # the total area

    # make sure TFR does not exceed y_upper
    TFR_modified = c()

    for( i in seq(1, length(TFR))) {
        if( TFR[i] > y_upper){
            TFR_modified[i] <- y_upper
        } else {
            TFR_modified[i] <- TFR[i]
        }
    }

    avrec = sum(diff(recall)*rollmean(log10(TFR_modified),2)) / sum_area

    if(is.na(avrec)){
        avrec = 0
    }

    avrec = round(avrec, digits=3)

    if(plots){

        unicolor = "darkblue"
        cex_main_size = 3.0

        #unicolor = "black"
        #cex_main_size = 3.5

        if(rerank){
            mainname = paste0("Motif Performance, AvRec=", avrec)
        } else{
            mainname = paste0("Dataset Performance, AvRec=", avrec)
            #mainname = paste0("Average Recall Curve, AvRec=", avrec)
        }

        # make sure that the upper border does not exceed y_upper
        mask_upper = TFR <= y_upper
        TFR_maskup = TFR[mask_upper]
        recall_maskup = recall[mask_upper]
        # plot the line when positives:negatives = 1:1
        png( filename = paste0(picname, ".png"), width = 800, height = 800 )
        par(oma=c(0,0,0,0), mar=c(6,6.5,5,1))
        plot(recall_maskup, TFR_maskup,
            main=mainname,
            log="y",
            xlim=c(0,1),ylim=c(y_lower,y_upper),
            xlab="", ylab="",
            type="l", lwd=7.5,
            col=unicolor,
            axes = FALSE, cex.main = cex_main_size
        )
        # color the area under the curve
        polygon(
            c(0, recall, 1),
            c(y_lower, pmin(TFR, y_upper), y_lower),
            col = convertcolor(unicolor,30),
            border = NA
        )

        mtext("Recall = TP/(TP+FN)", side=1, line=4.5, cex = cex_main_size)
        mtext("TP/FP Ratio", side=2, line=4, cex = cex_main_size)
        axis(1, at=c(0,0.5,1), labels = c(0,0.5,1), tick =TRUE, cex.axis=2.5, line=0)
        axis(2, at=c(y_lower,10,y_upper), labels = expression(10^0, 10^1, 10^2), tick=TRUE, cex.axis=2.5, line=0, las=1)
        text(x = 0.05,y = min(max(TFR)*1.1, 90), cex = 2.0, locator(), labels = c("1:1"), col=unicolor)

        if(max(TFR)>10){
            # calculate for the case when positives: negatives = 1:10
            # make sure that the lower border does not exceed y_lower
            mask_lower = TFR <= y_lower*10
            TFR_low = TFR[!mask_lower]
            recall_low = recall[!mask_lower]
            mask_2upper = TFR_low <= y_upper*10
            TFR_low = TFR_low[mask_2upper]
            recall_low = recall_low[mask_2upper]
            # plot the line when positives:negatives = 1:10
            par(new=T)
            plot(recall_low, TFR_low/10,
            main=mainname,
            log="y",
            xlim=c(0,1),ylim=c(y_lower,y_upper),
            xlab="", ylab="",
            type="l", lty=2, lwd=7.5,
            col=convertcolor(unicolor,50),
            axes = FALSE, cex.main = cex_main_size
            )
            # add label text
            text(x = 0.05, y = min(max(TFR/10)*1.1, 70), cex = 2.0, locator(), labels = c("1:10"), col=unicolor)
        }

        if(max(TFR)>100){
            # calculate for the case when positives: negatives = 1:100
            # make sure that the lower border does not exceed y_lower
            mask_lower2 = TFR <= y_lower*100
            TFR_low2 = TFR[!mask_lower2]
            recall_low2 = recall[!mask_lower2]
            mask_2upper2 = TFR_low2 <= y_upper*100
            TFR_low2 = TFR_low2[mask_2upper2]
            recall_low2 = recall_low2[mask_2upper2]

            # plot the line when positives:negatives = 1:100
            par(new=T)
            plot(recall_low2, TFR_low2/100,
            main=mainname,
            log="y",
            xlim=c(0,1),ylim=c(y_lower,y_upper),
            xlab="", ylab="",
            type="l", lty=3, lwd=7.5,
            col=convertcolor(unicolor,50),
            axes = FALSE, cex.main = cex_main_size
            )
            # add labels
            text(x = 0.05, y = min(max(TFR/100)*1.1, 55), cex = 2.0, locator(), labels = c("1:100"), col=unicolor)
        }

        box(lwd=2.5)
        invisible(dev.off())

        # export plots to .pdf file
        pdf( file = paste0(picname, ".pdf"), width = 10, height = 10 )
        par(oma=c(0,0,0,0), mar=c(6,6.5,5,1))
        plot(recall_maskup, TFR_maskup,
        main=mainname,
        log="y",
        xlim=c(0,1),ylim=c(y_lower,y_upper),
        xlab="", ylab="",
        type="l", lwd=7.5,
        col=unicolor,
        axes = FALSE, cex.main = cex_main_size
        )
        # color the area under the curve
        polygon(
        c(0, recall, 1),
        c(y_lower, pmin(TFR, y_upper), y_lower),
        col = convertcolor(unicolor,30),
        border = NA
        )

        mtext("Recall = TP/(TP+FN)", side=1, line=4.5, cex = 3.0)
        mtext("TP/FP Ratio", side=2, line=4, cex = cex_main_size)
        axis(1, at=c(0,0.5,1), labels = c(0,0.5,1), tick =TRUE, cex.axis=2.5, line=0)
        axis(2, at=c(y_lower,10,y_upper), labels = expression(10^0, 10^1, 10^2), tick=TRUE, cex.axis=2.5, line=0, las=1)
        text(x = 0.05,y = min(max(TFR)*1.1, 90), cex = 2.0, locator(), labels = c("1:1"), col=unicolor)

        if(max(TFR)>10){

            # plot the line when positives:negatives = 1:10
            par(new=T)
            plot(recall_low, TFR_low/10,
            main=mainname,
            log="y",
            xlim=c(0,1),ylim=c(y_lower,y_upper),
            xlab="", ylab="",
            type="l", lty=2, lwd=7.5,
            col=convertcolor(unicolor,50),
            axes = FALSE, cex.main = cex_main_size
            )
            # add label text
            text(x = 0.05, y = min(max(TFR/10)*1.1, 70), cex = 2.0, locator(), labels = c("1:10"), col=unicolor)
        }

        if(max(TFR)>100){
            # plot the line when positives:negatives = 1:100
            par(new=T)
            plot(recall_low2, TFR_low2/100,
            main=mainname,
            log="y",
            xlim=c(0,1),ylim=c(y_lower,y_upper),
            xlab="", ylab="",
            type="l", lty=3, lwd=7.5,
            col=convertcolor(unicolor,50),
            axes = FALSE, cex.main = cex_main_size
            )
            # add labels
            text(x = 0.05, y = min(max(TFR/100)*1.1, 55), cex = 2.0, locator(), labels = c("1:100"), col=unicolor)
        }

        box(lwd=2.5)
        invisible(dev.off())
    }
    # access the avrec value
    return( list(avrec=avrec) )
}

#--------------------------
#
# motif evaluation function
#
#--------------------------

evaluateMotif = function( pvalues, filename, rerank, data_eta0 ){

    if( rerank ){
        eta0set = NULL
        pn_pval <- paste0( filename, '_motifPval' )
        pn_rrc  <- paste0( filename, '_motifRRC' )
    } else {
        eta0set = data_eta0
        pn_pval <- paste0( filename, '_dataPval' )
        pn_rrc  <- paste0( filename, '_dataRRC' )
    }

    if(max(pvalues) > 1 | min(pvalues) <= 0){
        stop("Error: input p-values must all be in the range 0 to 1!")
    }

    # use fdrtool to estimate eta0
    result_fdrtool = fdrtool(pvalues, statistic="pvalue", plot=FALSE, eta0set=eta0set)

    # get the global fdr values and estimate of the weight eta0 of the null component
    fdr	    <- result_fdrtool$qval
    eta0 	<- result_fdrtool$param[3]

    if( eta0 >= 1 ){
        #stop("estimated eta0 >= 1. No positives in the input set.")
        eta0 = 0.9999
    }

    # plot p-value density plot
    if(plots) plotPvalStat(pvalues, filename=pn_pval, eta0=eta0, data_eta0=data_eta0, rerank)

    # calculate recall
    len         = length(fdr)
    list        <- seq(1, len)

    recall  <- ( 1 - fdr ) * list / ( 1 - eta0 ) / len

    # cut values when recall is larger than 1
    cutoff = len
    for(i in list){
        if( recall[i] > 1 ){
            cutoff = i-1
            break
        }
    }

    recall=recall[1:cutoff]
    fdr=fdr[1:cutoff]

    # set FP / TP ratio to 1:1
    mfold  = eta0 / (1-eta0)
    tfr <- numeric(len)

    if(max(fdr)>0){
        tfr <- (1-fdr)/fdr * mfold
    } else {
        tfr = 1
    }

    recall[cutoff] = 1
    tfr[cutoff] = 1

    if(eta0 < data_eta0) eta0 = data_eta0

    # plot TP/FP vs. recall plot
    rrc     = plotRRC(pn_rrc, recall, tfr, rerank)
    avrec   = rrc$avrec     # get AURRC score

    # return to the results
    return( list(avrec=avrec, eta0=eta0) )

}

#--------------
#
# main function
#
#--------------
file_suffix = ".zoops.stats"

results = c()
resultTitle = paste0(c("TF", "#", "d_avrec", "d_occur", "m_avrec", "m_occur"), collapse="\t")
results = c(results, resultTitle)

if( length(Sys.glob(paste(c(dir, "/", prefix, "*", file_suffix), collapse=""))) ==0 ){
    stop("Error: no input file exists in the folder!")
}

for (f in Sys.glob(paste(c(dir, "/", prefix, "*", file_suffix), collapse=""))) {

    front_end=paste0(c(dir, "/", prefix), collapse="")
    # get motif number from the filename; Note: important for motif reranking
    motif_id <- sub(front_end, "", f, fixed = TRUE)     # note: fixed=TRUE is very very important for recognizing special characters
    motif_id <- sub(file_suffix, "", motif_id, fixed = TRUE)
    motif_num <- unlist(strsplit(motif_id, "_motif_"))[-1]


    if( length(motif_num) == 0){
        motif_num = "NaN"
    }

    # get a filename for each motif
    filename = paste(c(dir, "/", prefix, motif_id), collapse="")

    # read in p-values from file
    first_row   <- read.table(f, nrows=1)
    stats       <- read.table(f, skip=1)
    pvalues     <- stats$V5
    #pvalues     <- unique(stats$V5)
    mfold       <- as.numeric(first_row[6])
    occurrence  <- as.numeric(first_row[7])

    # check if any p-values are missing
    if( sum(is.na(pvalues)) > 0 ){
        stop("Error: Some of the input p-values are missing!")
    }

    if( is.na(mfold) ){
        stop("Error: mfold value is missing in the input file!")
    }

    if( is.na(occurrence) ){
        stop("Error: Occurrence value is missing in the input file!")
    }

    # avoid the rounding errors when p-value = 0 or p-value > 1
    for(i in seq(1, length(pvalues))){
        if( pvalues[i] > 1 ){
            pvalues[i] = 1
        } else if (pvalues[i] == 0 ){
            pvalues[i] = 1e-10
        }
    }

    # calculate eta0 based on negative/postitive ratio
    data_eta0 = mfold / ( 1+mfold )

    # evaluate motif on the dataset
    eval_dataset    = evaluateMotif(pvalues, filename = filename, rerank=FALSE, data_eta0=data_eta0)
    data_occur      = round(occurrence, digits=4)               # acquire motif occurrence
    data_avrec      = eval_dataset$avrec                        # acquire AvRec score

    # evaluate motif indenpendent from dataset
    eval_motif      = evaluateMotif(pvalues, filename = filename, rerank=TRUE, data_eta0=data_eta0)
    motif_eta0      = eval_motif$eta0
    motif_occur     = round((1-motif_eta0)*(1+mfold), digits=4) # calculate motif occurrence
    motif_avrec     = eval_motif$avrec                          # acquire AvRec score

    # output the result
    resultString = paste0(c(prefix, motif_num, data_avrec, data_occur, motif_avrec, motif_occur), collapse="\t")

    if(web){
        message(motif_occur)
    }

    #print( resultString )

    results = c(results, resultString)
}
outConn <- file(ofile)
writeLines(results, outConn)
close(outConn)
