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
suppressMessages(library( argparse ))
# load "LSD" library for plotting
suppressMessages(library( LSD ))

suppressMessages(library( gdata ))

#-----------------------------
#
# Parameter setting
#
#-----------------------------
parser <- ArgumentParser(description="plot the motif distribution")

# positional arguments
parser$add_argument('maindir',  help="directory that contains the target file")
parser$add_argument('prefix',   help="prefix of the target file")

# optional arguments
parser$add_argument('--anchor', type="integer", help="Define where the starting point is")
parser$add_argument('--anchorname', type="character", default='0', help="Define the name of the starting point")
parser$add_argument('--ss', type="integer", help="Define whether motifs are searched on single strand or not")

args        <- parser$parse_args()
maindir     <- args$maindir
file_prefix <- args$prefix
anchor      <- args$anchor
anchorname  <- args$anchorname
ss          <- args$ss

#-----------------------------
#
# Read in the file
#
#-----------------------------
file_suffix = ".pi"

# plot parameters
png_width   <- 1000
png_height  <- 400
label_size  <- 2
lwd_size    <- 1.5

main_title  = "Motif Positional Prior"

# strip away the trailing slash if used by the user
maindir <- gsub('/$', '', maindir)

file_glob = paste(c(file_prefix, '*', file_suffix), collapse="")
full_glob = file.path(maindir, file_glob)
if( length(Sys.glob(full_glob)) == 0 ){
    stop("Error: no input file exists in the folder!")
}

for( file in Sys.glob(full_glob) ){

    # get motif number from the filename
    prefix   <- file.path(maindir, file_prefix)
    motif_id <- sub(prefix, "", file)
    motif_id <- sub(file_suffix, "", motif_id)

    # get a filename for each motif
    file_name   = paste(c(file_prefix, motif_id, file_suffix), collapse="")
    filename    = file.path(maindir, file_name)
    line_number = as.integer(system2("wc", args=c("-l", filename, " | awk '{print $1}'" ), stdout = TRUE))

    # read in the data
    table <- try(read.table(filename,
                            fileEncoding="latin1", as.is=TRUE, na.strings = "NA",
                            fill = TRUE, strip.white = TRUE, skip=1, sep = '\t'))

    motif_positions     = c(table$V1)
    whole_region        = line_number
    if( is.null(ss) ){
        strand_length   = whole_region / 2
        pos_positions   = motif_positions[0:strand_length]
        neg_positions   = tail(motif_positions, strand_length)

        # turn negative strand upside down
        neg_positions   <- -neg_positions
        neg_positions   <- rev(neg_positions)

        # set y-scale
        yscale          = c(min(neg_positions), max(pos_positions))

    } else {
        strand_length   = whole_region
        pos_positions   = motif_positions
        yscale          = c()
    }
    # define the starting point of the sequences
    if( is.null(anchor) ){
        strand_center   = strand_length / 2
    } else {
        strand_center   = anchor
    }

    interval_left   = strand_center
    interval_right  = strand_length - strand_center

    picname <- paste(c(file_prefix, motif_id, "_PositionalPrior.png"), collapse="")
    picname <- file.path(maindir, picname)
    png(filename = picname, width=png_width, height=png_height)
    par(oma=c(0,0,0,0), mar=c(6,6.5,5,2))

    plot(pos_positions,
        main=main_title,
        xlab="", ylab="",
        type="l",
        lwd=lwd_size*3,
        col="darkblue",
        axes=FALSE,
        cex.axis=label_size,
        cex.main=label_size,
        xlim = c(0, strand_length),
        ylim = yscale
    )
    legend("topright", legend="+ strand",
            col="darkblue", cex=label_size,
            bty="n", text.col="darkblue")
    if( is.null(ss) ){
        # plot positional priors on the negative strand
        lines(neg_positions,
                #type="l",
                lwd=lwd_size*3,
                col="darkred")
        legend("bottomright", legend="- strand",
            col="darkred", cex=label_size,
        bty="n", text.col="darkred")

    }

    mtext("Position relative to sequence center", side=1, line=4.5, cex=label_size)
    mtext("Positional prior", side=2, line=4, cex=label_size)
    abline(h=0, v=strand_center, col="grey", lwd=lwd_size*2)
    #axis(1, at=c(0, strand_center, strand_length),
    #labels = c(-interval_left, anchorname,  interval_right),
    axis(1, at=c(0, strand_center-10, strand_center, strand_center+10, strand_length),
    labels = c(-interval_left, -10, anchorname, 10, interval_right),
    tick = FALSE, cex.axis=label_size, line=1)
    axis(2, tick = FALSE, cex.axis=label_size, line=0.5)
    box(lwd=lwd_size)

    invisible(dev.off())

}
