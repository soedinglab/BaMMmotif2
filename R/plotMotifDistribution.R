#!/usr/bin/env Rscript

# Title     : plot the positional distribution of the query motif on the sequences
# Objective : plot the positional distribution of the query motif on the sequences
# Created by: Anja Kiesel and modified by Wanwan Ge
# Created on: 01.09.17


#-----------------------------
#
# R sources
#
#-----------------------------
# load "argparse" library for parsing arguments
library( argparse )
# load "LSD" library for plotting
library( LSD )

suppressMessages( library( gdata ) )

#-----------------------------
#
# Parameter setting
#
#-----------------------------
parser <- ArgumentParser(description="plot the motif distribution")
# positional arguments
parser$add_argument('target_directory', help="directory that contains the target file")
parser$add_argument('prefix', help="prefix of the target file")
# optional arguments
parser$add_argument("--ss", type="logical", default=FALSE, help="flag for searching only on single strand" )

args <- parser$parse_args()
maindir <- args$target_directory
basename <- args$prefix
revComp <- !args$ss

#-----------------------------
#
# Read in the file
#
#-----------------------------
positions <- read.table( paste0( maindir, '/', basename, ".positions" ),
                        fileEncoding="latin1", as.is=TRUE, na.strings = "NA",
                        fill = TRUE, strip.white = TRUE, skip=1)
all_positions = c()
pattern_range = c()
pattern_range = c( pattern_range, positions$V4 )
#pattern_range = c( pattern_range, positions$V7 )
seq.length = as.numeric( positions$V2[1] )
for( i in seq(1, length(pattern_range)) ){
    x <- unlist( regmatches(pattern_range[i], gregexpr('\\(?[0-9,]+', pattern_range[i])))
    x <- as.numeric(gsub('\\(', '-', gsub(',', '', x)))
    motif_pos = as.integer( ( x[2]+x[1] ) /2 )
    all_positions = c(all_positions, na.omit(motif_pos))
}

picname <- paste0( maindir, basename,"_distribution.jpeg")
jpeg( filename = picname, width = 800, height = 800, quality = 100 )
par(oma=c(0,0,0,0), mar=c(6,6.5,5,2))

if( revComp ){
    ################################################
    #   plot motif distribution on double strands  #
    ################################################
    pos.strand = density(all_positions[which(all_positions < seq.length)] )
    neg.strand = density(all_positions[which(all_positions >= seq.length)] - seq.length )
    # turn neg strand upside down
    neg.strand$y <- neg.strand$y*-1

    plot(pos.strand,
        main="Motif Positions",
        xlab="", ylab="",
        type="l", lwd=7.5,
        col="darkblue",
        axes=FALSE, cex.axis=3.0 ,3.0, cex.main=3.0,
        xlim = c(0,seq.length),
        ylim=c(min(neg.strand$y),max(pos.strand$y))

    )
    mtext("Position on Sequence", side=1, line=4.5, cex = 3.5)
    mtext("Density", side=2, line=4, cex = 3.5)
    axis(1, labels=FALSE)
    axis(1, tick = FALSE, cex.axis=3.0, line=1)
    axis(2, labels=FALSE)
    axis(2, tick = FALSE, cex.axis=3.0, line=0.5)
    polygon(pos.strand, col=convertcolor("darkblue",30),border = NA)
    lines(neg.strand, type="l", lwd=7.5, col="darkred")
    polygon(neg.strand, col=convertcolor("darkred",30),border = NA)
    legend("topright",legend = "+ strand",col = "darkblue",cex = 2.5, bty = "n",text.col = "darkblue")
    legend("bottomright",legend = "- strand",col = "darkred",cex = 2.5, bty = "n", text.col="darkred")
    abline(h=0)
    box(lwd=2.5)
    invisible(dev.off())

} else {
    ################################################
    #   plot motif distribution on single strand   #
    ################################################
    pos <- density(all_positions)
    plot(pos,
        main="Motif Positions",
        xlab="", ylab="",
        type="l", lwd=7.5,
        col="darkblue",
        axes=FALSE, cex.axis=3.0 ,3.0, cex.main=3.0
    )
    mtext("Position on Sequence", side=1, line=4.5, cex = 3.5)
    mtext("Density", side=2, line=4, cex = 3.5)
    axis(1, labels=FALSE)
    axis(1, tick = FALSE, cex.axis=3.0, line=1)
    axis(2, labels=FALSE)
    axis(2, tick = FALSE, cex.axis=3.0, line=0.5)
    polygon(pos, col=convertcolor("darkblue",30),border = NA)
    box(lwd=2.5)
    invisible(dev.off())
}