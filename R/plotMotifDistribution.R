#!/usr/bin/env Rscript

# Title     : plot the positional distribution of the query motif on the sequences
# Objective : plot the positional distribution of the query motif on the sequences
# Created by: Anja Kiesel and modified by Wanwan Ge
# Created on: 01.09.17

library( LSD )
suppressMessages( library( gdata ) )
################################################
#             read in parameters               #
################################################
#dir = "/home/kiesel/Desktop/BaMM_webserver/media/d55c520a-a9ff-4d2d-830e-4f811fc62fb8/"
#FileName = "positiveSequences_motif_1"
revComp= TRUE
#FastaFileName="positiveSequences"

args <- commandArgs(trailingOnly = TRUE)
splits <- strsplit(args,split='=')
for(split in splits) {
    if(split[1] == '--output_dir') {maindir <- split[2]}
    else if(split[1] == '--file_name_in'){FileName <- split[2]}
    else if(split[1] == '--revComp'){revComp <- split[2]}
    else if(split[1] == '--fasta_file_name'){FastaFileName <- split[2]}
}

plot_PRC = FALSE
plot_Dist = TRUE
if( plot_PRC ){
################################################
#                read in files                 #
################################################
infile_stats <- read.table( paste0( maindir, FileName, ".zoops.stats" ), fileEncoding="latin1", header = FALSE, skip=1)
colnames(infile_stats) <- c("TP","FP","FDR","Recall","p-value_20")
Recall <- infile_stats$Recall
Precision <- (1-infile_stats$FDR)
N = dim(infile_stats)[1]

################################################
#                  calculate AUC               #
################################################

auc <- sum(Precision)/N
message(round(auc,digits=3))

################################################
#                  plot curve                  #
################################################
picname <- paste0( maindir, FileName,"_FDR.jpeg")
jpeg( filename = picname, width = 800, height = 800, quality = 100 )
par(oma=c(0,0,0,0), mar=c(6,6.5,5,1))

plot(Recall, Precision,
    main="Motif Performance",
    xlim=range(0,1), ylim=range(0,1),
    xlab="", ylab="",
    type="l", lwd=7.5,
    col="darkblue",
    axes = FALSE, cex.main=3.0)
    mtext("Recall", side=1, line=4.5, cex = 3.5)
    mtext("Precision", side=2, line=4, cex = 3.5)
    axis(1, at=c(0,0.5,1),labels = c(0,0.5,1),tick = FALSE, cex.axis=3.0, line=1)
    axis(2, at=c(0,0.5,1),labels = c(0,0.5,1),tick = FALSE, cex.axis=3.0)
    polygon(c(0,Recall,1), c(0,Precision,0), col = convertcolor("darkblue",30),border = NA)
    text(x = 0.5,y = 0.3,labels = paste0("AUC = ",round(auc,digits=3)),cex = 3.5)
    box(lwd=2.5)
    invisible(dev.off())

    picname <- paste0( maindir, FileName,"_Positions.jpeg")
    jpeg( filename = picname, width = 800, height = 800, quality = 100 )
    par(oma=c(0,0,0,0), mar=c(6,6.5,5,2))
}

if( plot_Dist ){

################################################
#                read in files                 #
################################################

    position_file <- read.table( paste0( maindir, FileName, ".positions" ),
                                fileEncoding="latin1", header=TRUE )

    print(position_file$header)

    colnames(position_file) <- c(position_file$header)

    all_positions = c()
    for( i in c(2:dim(position_file)[2])){
        all_positions = c(all_positions, na.omit(position_file[,i]))
    }

if(revComp == "True"){
    # find out sequence length
    data = scan(paste0( maindir, FastaFileName ), what = list(fasta=""),quiet=TRUE)
    next_header = grep( x = data[[1]], pattern=">")[2]
    seq.sample = paste0(data[[1]][2:(next_header-1)],collapse="")
    seq.length = nchar(seq.sample)

    print(seq.length)

    pos.strand = density(all_positions[which(all_positions < seq.length)])
    neg.strand = density(all_positions[which(all_positions >= seq.length)]-seq.length)


    #    pos.strand = density(c(na.omit(positions[which(positions[,2] < seq.length),2]),na.omit(positions[which(positions[,3] < seq.length),3])))
    #    neg.strand = density(c(na.omit(positions[which(positions[,3] >= seq.length),3]-seq.length),na.omit(positions[which(positions[,3] >= seq.length),3]-seq.length)))

    # turn neg strand upside down
    neg.strand$y <- neg.strand$y*-1

    plot(pos.strand,
        main="Motif Positions",
        xlab="", ylab="",
        type="l", lwd=7.5,
        col="darkblue",
        axes=FALSE, cex.axis=3.0 ,3.0, cex.main=3.0,
        xlim = c(0,seq.length), ylim=c(min(neg.strand$y),max(pos.strand$y)))
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
    pos <- density(all_positions)
    #    pos <- density(c(na.omit(positions[,2]),na.omit(positions[,3])))
    ################################################
    #                  plot curve                  #
    ################################################

    plot(pos,
    main="Motif Positions",
    xlab="", ylab="",
    type="l", lwd=7.5,
    col="darkblue",
    axes=FALSE, cex.axis=3.0 ,3.0, cex.main=3.0)
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


################################################
#            calculate Occurrence              #
################################################

occurrence = sum(!is.na(position_file[,2]))/dim(position_file)[1]
message(round(occurrence,digits=3))
}