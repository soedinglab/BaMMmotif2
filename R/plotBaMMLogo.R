#!/usr/bin/env Rscript
#...............................................................................
# 
# Documentations
#
#...............................................................................

### This R script is used for plotting PWMs and higher-order sequence logo from 
### bamm-formatted files.
### It was written by Dr. Matthias Siebert, 
### and modified by Wanwan Ge (wanwan.ge@mpibpc.mpg.de)
### copy right 2017 Soeding Group

### Input file
# * Plot higher-order sequence logo from bamm-formatted file.
# e.g. basename.ihbcp (conditional probability file)
# or   basename.ihbp  (frequence file)

### Result
# * Higher-order sequence logo
# e.g. basename.pdf

### command line
# plotBaMMlogo.R path_to_input_bamm_file logo_order
# Note: The logo order can be from 0 to 5
# e.g.
# plotBaMMlogo.R /tmp/exampleset/example.ihbcp 0 

#...............................................................................
#
# Libraries
#
#...............................................................................

# A package for parsing arguments
library( argparse )

#...............................................................................
#
# Function: read in BaMM file
#
#...............................................................................

readBaMM <- function( bamm_filename ){
  
  list_flatprobs <- readLines( bamm_filename )  # character vector
  list_flatprobs <- lapply( strsplit( list_flatprobs, " " ), as.double ) # double list
  
  # separate model positions
  blanks <- which( sapply( list_flatprobs, length ) == 0 )
  
  # check if same order at single positions
  if( length( unique( blanks[-1] - blanks[-length( blanks )] ) ) > 1 ){
    stop( paste( "Incorrect IMM flat file format: ", bamm_filename, sep="" ) )
  }
  
  # remove leading blank (if any)
  if( blanks[1] == 1 ){
    list_flatprobs <- list_flatprobs[-1]
    blanks <- blanks[-1]
  }
  
  l <- length( list_flatprobs )
  # remove trailing blank (if any)
  if( rev( blanks )[1] == l ){
    list_flatprobs <- rev( rev( list_flatprobs )[-1] )
    blanks <- rev( rev( blanks )[-1] )
  }
  
  # calculate model positions
  L <- length( which( sapply( list_flatprobs, length ) == 0 ) ) + 1
  # calculate max. model order
  maxorder <- ( ( length( list_flatprobs )-( L-1 ) ) / L ) - 1
  
  MAXORDER <- 9 # as implemented so far
  if( maxorder > MAXORDER ){
    stop( paste( "Max. order: ", MAXORDER, sep="" ) )
  }
  
  list_probs <- list() # initialize IMM data structure
  if( maxorder >= 0 ) list_probs[[ 1]] <- array( NA, c( L, rep( 4,  1 ) ) )
  if( maxorder >= 1 ) list_probs[[ 2]] <- array( NA, c( L, rep( 4,  2 ) ) )
  if( maxorder >= 2 ) list_probs[[ 3]] <- array( NA, c( L, rep( 4,  3 ) ) )
  if( maxorder >= 3 ) list_probs[[ 4]] <- array( NA, c( L, rep( 4,  4 ) ) )
  if( maxorder >= 4 ) list_probs[[ 5]] <- array( NA, c( L, rep( 4,  5 ) ) )
  if( maxorder >= 5 ) list_probs[[ 6]] <- array( NA, c( L, rep( 4,  6 ) ) )
  if( maxorder >= 6 ) list_probs[[ 7]] <- array( NA, c( L, rep( 4,  7 ) ) )
  if( maxorder >= 7 ) list_probs[[ 8]] <- array( NA, c( L, rep( 4,  8 ) ) )
  if( maxorder >= 8 ) list_probs[[ 9]] <- array( NA, c( L, rep( 4,  9 ) ) )
  if( maxorder >= 9 ) list_probs[[10]] <- array( NA, c( L, rep( 4, 10 ) ) )
  
  i <- 1     # flatprobs entries
  pos <- 1   # model positions
  order <- 0 # model orders
  while( i <= l ){
    
    if( order > maxorder ){ # blank (next model position)
      pos <- pos + 1
      order <- 0
    } else {                # no blank (next order at current model position)
      switch( order + 1,
              list_probs[[ 1]][pos,] <- list_flatprobs[[i]],
              list_probs[[ 2]][pos,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  2 ) ),  2:1 ),
              list_probs[[ 3]][pos,,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  3 ) ),  3:1 ),
              list_probs[[ 4]][pos,,,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  4 ) ),  4:1 ),
              list_probs[[ 5]][pos,,,,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  5 ) ),  5:1 ),
              list_probs[[ 6]][pos,,,,,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  6 ) ),  6:1 ),
              list_probs[[ 7]][pos,,,,,,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  7 ) ),  7:1 ),
              list_probs[[ 8]][pos,,,,,,,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  8 ) ),  8:1 ),
              list_probs[[ 9]][pos,,,,,,,,,] <- aperm( array( list_flatprobs[[i]], rep( 4,  9 ) ),  9:1 ),
              list_probs[[10]][pos,,,,,,,,,,] <- aperm( array( list_flatprobs[[i]], rep( 4, 10 ) ), 10:1 ),
      )
      order <- order + 1;
    }
    i <- i + 1
    
  }
  
  return( list_probs )
}

#...............................................................................
#
# Source
#
#...............................................................................
source( "generateLogo.R" )

#...............................................................................
#
# Parameter setting
#
#...............................................................................

parser <- ArgumentParser(description='plot higher-order bamm sequence logo')
parser$add_argument('input_path', help='input path to the bamm file')
parser$add_argument('logo_order', help='logo order, it must be not higher than bamm model order')
args <- parser$parse_args()

input_path <- args$input_path

# directory to hold the output log file
dirname <- dirname(input_path)

# read in filename of the bamm files without extension
# probabilities (.ihbp) mandatory
# conditional probabilities (.ihbcp) and background frequencies (.hbp) optional
basename <- basename(input_path)
splits <- strsplit( basename, split = '\\.' )[[1]]
length <- length(splits)
filename <- splits[1]
for( i in seq(2, length-1) ){
  filename <- paste( filename, splits[i], sep='.')
}
# read in logo order (max. order: 5)
order <- args$logo_order

# background sequence logo
bg <- FALSE

# use background (instead of uniformely distributed) nucleotide frequencies for
# 0th-order 1-mer contributions (unapplicable for background sequence logos)
useFreqs <- FALSE

# logarithm base to use in calculating k-mer contributions
base <- 2

# information content (instead of probability) column scaling
icColumnScale <- TRUE
# information content (instead of probability) letter scaling
icLetterScale <- TRUE

# further logo display parameters
xaxis <- TRUE
yaxis <- TRUE
xfontsize <- 15
yfontsize <- 15
alpha <- .3 # icLetterScale only

# .pdf parameters
width <- 7


#...............................................................................
#
# Plot higher-order sequence logo(s)
#
#...............................................................................

if( bg ){
    png( file.path( dirname, paste( filename, "-bg-logo-order-", paste( order,
         collapse="-" ), ifelse( icColumnScale, "-icColumnScale", "" ), ifelse(
         icLetterScale, "-icLetterScale", "" ), ".png", sep="" ) ) )
    hoBgSeqLogo( filename=file.path( dirname, filename ), order=order, base=base,
                 icColumnScale=icColumnScale, icLetterScale=icLetterScale,
                 xaxis=xaxis, yaxis=yaxis, xfontsize=xfontsize,
                 yfontsize=yfontsize, alpha=alpha )
} else {
    png( file.path( dirname, paste( filename, "-logo-order-", order, ".png", sep="" ) ) )
    hoSeqLogo( filename=file.path( dirname, filename ), order=order,
               useFreqs=useFreqs, base=base, icColumnScale=icColumnScale,
               icLetterScale=icLetterScale, xaxis=xaxis, yaxis=yaxis,
               xfontsize=xfontsize, yfontsize=yfontsize, alpha=alpha )
}
dev.off()
