#...................................................................................................
#
# Documentations
#...................................................................................................

# Job
# * Read in (inhomogeneous) IMM saved as flat file
#
# Results
# * The (inhomogeneous) IMM (list of 1(+1) to max.order+1(+1) dimensional arrays)


#...................................................................................................
#
# Parameters
#...................................................................................................

# file_IMM: (inhomogeneous) IMM flat file


#...................................................................................................
#
# Functions
#...................................................................................................

Utils.readIMMs.flat <- function( file_IMM ){

    list_flatprobs <- readLines( file_IMM ) # character vector
    list_flatprobs <- lapply( strsplit( list_flatprobs, " " ), as.double ) # double list

    # separate model positions
    blanks <- which( sapply( list_flatprobs, length ) == 0 )
    # check if same order at single positions
    if( length( unique( blanks[-1] - blanks[-length( blanks )] ) ) > 1 ){
	stop( paste( "Incorrect IMM flat file format: ", file_IMM, sep="" ) )
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

    i <- 1 # flatprobs entries
    pos <- 1 # model positions
    order <- 0 # model orders
    while( i <= l ){

	if( order > maxorder ){ # blank (next model position)
	    pos <- pos + 1
	    order <- 0
	}
	else{ # no blank (next order at current model position)
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
