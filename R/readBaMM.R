################################################################################
#
# Documentation
# 
################################################################################

# Read in (homogeneous) BaMM from file with BaMM probabilities as obtained from
# BaMM!motif

################################################################################
#
# Code to read in (homogeneous) BaMM from file
# 
################################################################################

readBaMM <- function( file ){

    list_flatprobs <- readLines( file ) # character vector
    list_flatprobs <- list_flatprobs[!( grepl( "^#", list_flatprobs ) )] # remove comment lines
    list_flatprobs <- lapply( strsplit( list_flatprobs, " " ), as.double ) # double list

    # separate model positions (using blank lines)
    blanks <- which( sapply( list_flatprobs, length ) == 0 )
    # check whether each position has the same order (i.e. number of lines)
    if( length( unique( blanks[-1] - blanks[-length( blanks )] ) ) > 1 ){
        stop( paste( "Wrong BaMM file format: ", file, sep="" ) )
    }

    if( length( blanks ) > 0 ){
        # remove leading blanks (if any)
        if( blanks[1] == 1 ){
            list_flatprobs <- list_flatprobs[-1]
            blanks <- blanks[-1]
        }
    }

    if( length( blanks ) > 0 ){
        # remove trailing blank (if any)
        if( rev( blanks )[1] == length( list_flatprobs ) ){
            list_flatprobs <- rev( rev( list_flatprobs )[-1] )
            blanks <- rev( rev( blanks )[-1] )
        }
    }

    # calculate model positions
    L <- length( which( sapply( list_flatprobs, length ) == 0 ) ) + 1
    # calculate maximum model order
    maxorder <- ( ( length( list_flatprobs )-( L-1 ) ) / L ) - 1

    MAXORDER <- 9 # as implemented so far
    if( maxorder > MAXORDER ){
        stop( paste( "Maximum BaMM order implemented so far: ", MAXORDER, sep=""
              ) )
    }

    list_probs <- list()
    # initialize BaMM
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

    pos <- 1 # model positions
    order <- 0 # model orders
    for( i in 1:length( list_flatprobs ) ){

        if( order > maxorder ){ # blank (next model position)
            pos <- pos + 1
            order <- 0
        } else{ # no blank (next order at current model position)
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
    }
    return( list_probs )
}
