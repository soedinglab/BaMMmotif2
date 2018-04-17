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

args        <- parser$parse_args()
maindir     <- args$maindir
file_prefix <- args$prefix

#-----------------------------
#
# Read in the file
#
#-----------------------------
file_suffix = ".occurrence"

if( length(Sys.glob(paste(c(maindir, '/', file_prefix, "*", file_suffix), collapse=""))) == 0 ){
    stop("no input file exists in the folder!")
}

for( file in Sys.glob(paste(c(maindir, '/', file_prefix, "*", file_suffix), collapse="")) ){

    # get motif number from the filename
    motif_id <- sub(paste(c(maindir, '/', file_prefix), collapse=""), "", file)
    motif_id <- sub(file_suffix, "", motif_id)

    # get a filename for each motif
    filename = paste0(c(maindir, file_prefix, motif_id, file_suffix), collapse="")
    line_number = as.integer(system2("wc", args=c("-l", filename, " | awk '{print $1}'" ), stdout = TRUE))

    label_size = 2
    main_title = "Motif Distribution"
    if( line_number < 3 ){
        print("The input file is empty. No query motif is found in the sequence set.")
        # print out an empty image
        picname <- paste0(maindir, file_prefix, motif_id, "_distribution.png")
        png(filename=picname, width=1000, height=400)
        par(oma=c(0,0,0,0), mar=c(6,6.5,5,2))
        plot(y=0,
            main=main_title,
            xlab="", ylab="",
            type="l", lwd=7.5,
            col="darkblue",
            axes=FALSE, cex.axis=label_size, cex.main=label_size+0.5,
            xlim = c(-100, 100)
        )
        abline(v=0, col="grey", lwd=3)
        abline(h=0, col="grey", lwd=3)
        mtext("Position relative to peak summit", side=1, line=4.5, cex=label_size)
        mtext("Density", side=2, line=4, cex=label_size)

        axis(1, at = c(-100, 0, 100),
            labels = c(-100, 0, 100),
            tick = FALSE, cex.axis=label_size, line=1)
            axis(2, tick = FALSE, cex.axis=label_size, line=0.5)

        legend("topright",legend="no motif found", col="darkblue", cex=label_size, bty="n", text.col="darkblue")
        box(lwd=2.5)
        invisible(dev.off())

    } else {
        # read in the data
        table <- try(read.table(filename,
                            fileEncoding="latin1", as.is=TRUE, na.strings = "NA",
                            fill = TRUE, strip.white = TRUE, skip=1, sep = '\t'))

        seqCount            = line_number - 1
        strand_length       = c(table$V2)
        max_strand_length   = max(strand_length)
        strand_center       = max_strand_length / 2
        strand_ind          = c(table$V3)
        negCount            = length(grep('-', strand_ind)) # count negative strands
        posCount            = seqCount - negCount           # count positive strands
        start               = as.numeric(lapply(strsplit(table$V4, split='..', fixed=TRUE), `[`, 1))
        end                 = as.numeric(lapply(strsplit(table$V4, split='..', fixed=TRUE), `[`, 2))
        motif_pos           = as.integer(( start+end ) / 2) - strand_length / 2 + strand_center

        pos_positions       <- rep(NA, posCount)
        neg_positions       <- rep(NA, negCount)

        if( negCount!= 0 ){
            pos_idx = 1
            neg_idx = 1
            for( i in seq(1, seqCount) ){
                if( strand_ind[i] == '+' ){
                    pos_positions[pos_idx] = motif_pos[i]
                    pos_idx = pos_idx + 1
                } else {
                    neg_positions[neg_idx] = motif_pos[i] - strand_length[i]
                    neg_idx = neg_idx + 1
                }
            }
        } else {
            pos_positions = motif_pos
        }

        interval = max(max(pos_positions)-strand_center, strand_center-min(pos_positions))
        interval = (as.integer(interval/10)+1) * 10
        picname <- paste0( maindir, file_prefix, motif_id, "_distribution.png")

        png(filename = picname, width = 800, height = 800)
        par(oma=c(0,0,0,0), mar=c(6,6.5,5,2))

        # calculate weights for kernel density estimation
        ## use 'counts / n' as weights:
        #pos.strand = density(unique(pos_positions), weights=table(pos_positions)/length(pos_positions))

        pos.strand = density((pos_positions))
        if(negCount!= 0){
            ## use 'counts / n' as weights:
            #neg.strand = density(unique(neg_positions), weights=table(neg_positions)/length(neg_positions))

            neg.strand = density((neg_positions))
            # turn neg strand upside down
            neg.strand$y <- neg.strand$y*-1

            # for plotting distribution on postive strand
            plot(pos.strand,
                main=main_title,
                xlab="", ylab="",
                type="l", lwd=7.5,
                col="darkblue",
                axes=FALSE, cex.axis=label_size, cex.main=label_size+0.5,
                xlim = c(strand_center - interval, strand_center + interval),
                ylim = c(min(neg.strand$y), max(pos.strand$y)
                )
            )

            # for plotting distribution on negative strand
            lines(neg.strand, type="l", lwd=7.5, col="darkred")
            polygon(neg.strand, col=convertcolor("darkred",30), border = NA)
            legend("bottomright",legend="- strand", col="darkred", cex=label_size, bty="n", text.col="darkred")

        } else {
            plot(pos.strand,
            main=main_title,
            xlab="", ylab="",
            type="l", lwd=7.5,
            col="darkblue",
            axes=FALSE, cex.axis=label_size, cex.main=label_size+0.5,
            xlim = c(strand_center - interval, strand_center + interval)
            )
        }

        abline(h=0, v=strand_center, col="grey", lwd=3)
        mtext("Position relative to peak summit", side=1, line=4.5, cex=label_size)
        mtext("Density", side=2, line=4, cex=label_size)
        axis(1, at=c(strand_center - interval, strand_center, strand_center + interval),
            labels = c(-interval, 0, interval),
            tick = FALSE, cex.axis=label_size, line=1)
        axis(2, tick = FALSE, cex.axis=label_size, line=0.5)
        polygon(pos.strand, col=convertcolor("darkblue", 30), border = NA)
        legend("topright",legend="+ strand", col="darkblue", cex=label_size, bty="n", text.col="darkblue")
        box(lwd=2.5)

        invisible(dev.off())
    }
}