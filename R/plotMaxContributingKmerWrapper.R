################################################################################
#
# Documentation
# 
################################################################################

# Plot the k-mer that contributes most to the information content
# using files with homogeneous BaMM probabilities as obtained from BaMM!motif

################################################################################
#
# R Sources
# 
################################################################################

# main directory
mainDir <- file.path( "", "home", "user", "workspace", "BaMMmotif-dev", "R" )

# function to plot max. contributing k-mer
source( file.path( mainDir, "plotMaxContributingKmer.R" ) )

################################################################################
#
# Parameter settings
# 
################################################################################

# directory with homogeneous BaMM files
bammDir <- file.path( mainDir, "BaMMs" )

# directory for resulting plot
plotDir <- file.path( mainDir, "Plots" )

# file with homogeneous BaMM probabilities as obtained from BaMM!motif
# omit filename extensions
filename <- "genome"

# homogeneous BaMM logo order (maximum order = 6)
order <- 2

# base of logarithm used to compute the information content
base <- 2

# use RNA alphabet (ACGU) instead of DNA (ACGT) alphabet
RNASeqLogo <- FALSE

# filename extension for homogeneous BaMM probabilities
file.extension.probs <- "hbp"

# filename extension for homogeneous BaMM conditional probabilities
file.extension.conds <- "hbcp"

# the width of the svg device in inches
width <- 7

# the height of the svg device in inches
height <- 7

# the default pointsize of plotted text (in big points)
pointsize <- 12

################################################################################
#
# Code to plot maximum contributing k-mer
# 
################################################################################

if( !( file.exists( plotDir ) ) ){
    dir.create( plotDir )
}

pdf( file.path( plotDir, paste( filename, "_max-contr-", order, "-mer.pdf",
                sep="" ) ), width=width, height=height, pointsize=pointsize )
opar <- par( no.readonly=TRUE )

plotMaxContributingBgKmer( filename=file.path( bammDir, filename ), order=order,
                           base=base, rna=RNASeqLogo,
                           probs.ending=file.extension.probs,
                           conds.ending=file.extension.conds )

par( opar )
dev.off()
