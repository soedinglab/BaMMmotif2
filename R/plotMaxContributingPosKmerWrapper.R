################################################################################
#
# Documentation
# 
################################################################################

# Plot positional k-mers that contribute most to the information content
# using files with BaMM probabilities as obtained from BaMM!motif

################################################################################
#
# R Sources
# 
################################################################################

# main directory
mainDir <- file.path( "", "home", "user", "workspace", "BaMMmotif-dev", "R" )

# function to plot max. contributing k-mers
source( file.path( mainDir, "plotMaxContributingKmer.R" ) )

################################################################################
#
# Parameter settings
# 
################################################################################

# directory with BaMM files
bammDir <- file.path( mainDir, "BaMMs" )

# directory for resulting plot
plotDir <- file.path( mainDir, "Plots" )

# file with BaMM probabilities as obtained from BaMM!motif
# omit filename extensions
filename <- "wgEncodeSydhTfbsHepg2Mafkab50322IggrabAlnRep0_summits102_restr5000_k2-K2-a1-20-60"

# BaMM logo order (maximum order = 5)
order <- 2

# base of logarithm used to compute the information content
base <- 2

# use RNA alphabet (ACGU) instead of DNA (ACGT) alphabet
RNASeqLogo <- FALSE

# use mononucleotide frequencies from <filename>.freqs file to calculate
# 0'th-order background frequencies (default: uniformly distributed nucleotide
# frequencies)
useFreqs <- FALSE

# filename extension for BaMM probabilities
file.extension.probs <- "probs"

# filename extension for BaMM conditional probabilities
file.extension.conds <- "conds"

# filename extension for BaMM mononucleotide background frequencies
file.extension.freqs <- "freqs"

# the width of the svg device in inches
width <- 7

# the height of the svg device in inches
height <- 7

# the default pointsize of plotted text (in big points)
pointsize <- 12

################################################################################
#
# Code to plot maximum contributing k-mer(s) of (homogeneous) BaMM
# 
################################################################################

if( !( file.exists( plotDir ) ) ){
    dir.create( plotDir )
}

pdf( file.path( plotDir, paste( filename, "_max-contr-", order, "-mer.pdf",
                sep="" ) ), width=width, height=height, pointsize=pointsize )
opar <- par( no.readonly=TRUE )

plotMaxContributingKmer( filename=file.path( bammDir, filename ), order=order,
                         base=base, rna=RNASeqLogo, useFreqs=useFreqs,
                         probs.ending=file.extension.probs,
                         conds.ending=file.extension.conds,
                         freqs.ending=file.extension.freqs )

par( opar )
dev.off()
