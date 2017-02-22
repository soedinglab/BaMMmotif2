################################################################################
#
# Documentation
# 
################################################################################

# Plot higher-order (homogeneous) BaMM logos
# using files with BaMM probabilities as obtained from BaMM!motif

# This code is built on code from the Bioconductor package seqLogo
# bioconductor.org/packages/release/bioc/html/seqLogo.html

################################################################################
#
# R Libraries
# 
################################################################################

library( grid )

################################################################################
#
# R Sources
# 
################################################################################

source( "readBaMM.R" )

################################################################################
#
# R functions to plot BaMM logo
# 
################################################################################

letterA <- function( x.pos, y.pos, ht, wt, id=NULL ){
	
    x <- c( 0,  4,  6, 10, 8, 5.0, 2, 0, 3.2, 3.6, 6.4, 6.8, 3.2 )
    y <- c( 0, 10, 10,  0, 0, 7.5, 0, 0, 3.0, 4.0, 4.0, 3.0, 3.0 )
    
    x <- 0.1 * x
    y <- 0.1 * y
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    if( is.null( id ) ){
        id <- c( rep( 1, 9 ), rep( 2, 4 ) )
    } else{
        id <- c( rep( id, 9 ), rep( id+1, 4 ) )
    }
    
    fill <- c( "#008000", "#008000" )
    col <- c( "#008000", "#008000" )
    
    list( x=x, y=y, id=id, fill=fill, col=col )
}

letterC <- function( x.pos, y.pos, ht, wt, id=NULL ){

    angle1 <- seq( 0.3+pi/2, pi, length=100 )
    angle2 <- seq( pi, pi*1.5, length=100 )
    
    x.l1 <- 0.5 + 0.5*sin( angle1 )
    y.l1 <- 0.5 + 0.5*cos( angle1 )
    x.l2 <- 0.5 + 0.5*sin( angle2 )
    y.l2 <- 0.5 + 0.5*cos( angle2 )

    x.l <- c( x.l1, x.l2 )
    y.l <- c( y.l1, y.l2 )
    
    x <- c( x.l, rev( x.l ) )
    y <- c( y.l, 1-rev( y.l ) )
    
    x.i1 <- 0.5 + 0.35*sin( angle1 )
    y.i1 <- 0.5 + 0.35*cos( angle1 )

    x.i1 <- x.i1[ y.i1 <= max( y.l1 ) ]
    y.i1 <- y.i1[ y.i1 <= max( y.l1 ) ]

    y.i1[1] <- max( y.l1 )
    
    x.i2 <- 0.5 + 0.35*sin( angle2 )
    y.i2 <- 0.5 + 0.35*cos( angle2 )
    
    x.i <- c( x.i1, x.i2 )
    y.i <- c( y.i1, y.i2 )
    
    x1 <- c( x.i, rev( x.i ) )
    y1 <- c( y.i, 1-rev( y.i ) )
    
    x <- c( x, rev( x1 ) )
    y <- c( y, rev( y1 ) )
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    if( is.null( id ) ){
        id <- rep( 1, length( x ) )
    } else{
        id <- rep( id, length( x ) )
    }
    
    fill <- "#0000FF"
    col <- "#0000FF"
    
    list( x=x, y=y, id=id, fill=fill, col=col )
}

letterG <- function( x.pos, y.pos, ht, wt, id=NULL ){

    angle1 <- seq( 0.3 + pi/2, pi, length=100 )
    angle2 <- seq( pi, 1.5*pi, length=100 )

    x.l1 <- 0.5 + 0.5*sin( angle1 )
    y.l1 <- 0.5 + 0.5*cos( angle1 )
    x.l2 <- 0.5 + 0.5*sin( angle2 )
    y.l2 <- 0.5 + 0.5*cos( angle2 )
    
    x.l <- c( x.l1, x.l2 )
    y.l <- c( y.l1, y.l2 )
    
    x <- c( x.l, rev( x.l ) )
    y <- c( y.l, 1-rev( y.l ) )
    
    x.i1 <- 0.5 + 0.35*sin( angle1 )
    y.i1 <- 0.5 + 0.35*cos( angle1 )

    x.i1 <- x.i1[ y.i1 <= max( y.l1 ) ]
    y.i1 <- y.i1[ y.i1 <= max( y.l1 ) ]

    y.i1[1] <- max( y.l1 )
    
    x.i2 <- 0.5 + 0.35*sin( angle2 )
    y.i2 <- 0.5 + 0.35*cos( angle2 )
    
    x.i <- c( x.i1, x.i2 )
    y.i <- c( y.i1, y.i2 )
    
    x1 <- c( x.i, rev( x.i ) )
    y1 <- c( y.i, 1-rev( y.i ) )
    
    x <- c( x, rev( x1 ) )
    y <- c( y, rev( y1 ) )
    
    h1 <- max( y.l1 )
    r1 <- max( x.l1 )
    
    h1 <- 0.4
    x.add <- c( r1, 0.5, 0.5, r1-0.2, r1-0.2, r1, r1 )
    y.add <- c( h1, h1, h1-0.1, h1-0.1, 0, 0, h1 )
    
    if( is.null( id ) ){
        id <- c( rep( 1, length( x ) ), rep( 2, length( x.add ) ) )
    } else{
        id <- c( rep( id, length( x ) ), rep( id+1, length( x.add ) ) )
    }
    
    x <- c( rev( x ), x.add )
    y <- c( rev( y ), y.add )
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    
    fill <- c( "#FFA500", "#FFA500" )
    col  <- c( "#FFA500", "#FFA500" )
    
    list( x=x, y=y, id=id, fill=fill, col=col )
}

letterT <- function( x.pos, y.pos, ht, wt, id=NULL ){
	
    x <- c( 0, 10, 10, 6, 6, 4, 4, 0 )
    y <- c( 10, 10, 9, 9, 0, 0, 9, 9 )

    x <- 0.1*x
    y <- 0.1*y
    
    x <- x.pos + wt*x
    y <- y.pos + ht*y
    
    if( is.null( id ) ){
        id <- rep( 1, 8 )
    } else{
        id <- rep( id, 8 )
    }
    
    fill <- "#FF0000"
    col <- "#FF0000"
    
    list( x=x, y=y, id=id, fill=fill, col=col )
}

letterU <- function( x.pos, y.pos, ht, wt, id=NULL ){

    w <- 2*pi/360
    i <- 0:360

    xInnerCircle <- ( cos( i*w )*3 )+5
    yInnerCircle <- ( sin( i*w )*3 )+5
    xInnerCircle <- c( 2, xInnerCircle[180:360], 8 )
    yInnerCircle <- c( 5, yInnerCircle[180:360], 5 )

    xOuterCircle <- ( cos( i*w )*5 )+5
    yOuterCircle <- ( sin( i*w )*5 )+5
    xOuterCircle <- c( 0, xOuterCircle[180:360], 10 )
    yOuterCircle <- c( 5, yOuterCircle[180:360], 5 )

    x <- c( 0, 2, 2, xInnerCircle, 8, 8, 10, rev( xOuterCircle ) )
    y <- c( 9, 9, 5, yInnerCircle, 5, 9, 9, rev( yOuterCircle ) )

    x <- 0.1*x
    y <- 0.1*y

    x <- x.pos + wt*x
    y <- y.pos + ht*y

    if( is.null( id ) ){
        id <- rep( 1, length( x ) )
    }else{
        id <- rep( id, length( x ) )
    }

    fill <- "#FF0000"
    col <-  "#FF0000"

    list( x=x, y=y, id=id, fill=fill, col=col )
} 

addLetter <- function( letters, which, x.pos, y.pos, ht, wt ){

    if( which == "A" ){ # which corresponds to letter
        letter <- letterA( x.pos, y.pos, ht, wt )
    } else if( which == "C" ){
        letter <- letterC( x.pos, y.pos, ht, wt )
    } else if( which == "G" ){
        letter <- letterG( x.pos, y.pos, ht, wt )
    } else if( which == "T" ){
        letter <- letterT( x.pos, y.pos, ht, wt )
    } else if( which == "U" ){
        letter <- letterU( x.pos, y.pos, ht, wt )
    } else{
        stop( "which must be A, C, G, T, or U" )
    }

    letters$x <- c( letters$x, letter$x )
    letters$y <- c( letters$y, letter$y )
    
    lastID <- ifelse( is.null( letters$id ), 0, max( letters$id ) )
    letters$id <- c( letters$id, lastID + letter$id )
    letters$fill <- c( letters$fill, letter$fill )
    letters$col <- c( letters$col, letter$col )

    letters
}

# calculate the information content for vector <x> using the logarithm to the
# base <base> (default: 2).
informationContent <- function( x, base=2 ){
    ifelse( all( x > 0 ), 2 + sum( x * log( x, base ) ), 0 )
}

# calculate the position-specific information content for the position-specific
# weight matrix <matrix> having positions in rows and nucleotide frequencies in
# columns (<margin>: 1) or vice versa (<margin>: 2). The information content is
# computed by using the logarithm to the base <base>.
informationContent.matrix <- function( matrix, margin=c( 1, 2 ), base=2 ){
    return( apply( matrix, margin[1], informationContent, base ) )
}

# calculate the higher-order position-specific information content of the
# inhomogeneous BaMM from k-mer contributions calculated by
# <kmerContributionsToInformationContent>. Discriminate between positive and
# negative contributions (<unsigned>: FALSE) or not (<unsigned>: TRUE).
informationContent.kmerContributions <- function( kmerContributions,
                                                  unsigned=TRUE ){

    maxorder <- length( kmerContributions )-1
    # n(ot) n(ull) indices
    nnindices <- which( !( sapply( kmerContributions, is.null ) ) )
    L <- unique( sapply( kmerContributions[ nnindices ], ncol ) )-1
    if( length( L ) > 1 ){
        stop( "Wrong kmerContributions format" )
    }

    if( unsigned ){
        ic <- array( NA, c( maxorder+1, L ) )
        
        for( i in 1:( maxorder+1 ) ){
            if( !( is.null( kmerContributions[[i]] ) ) ){
                ic[i,i:L] <- apply( kmerContributions[[i]][,-(1:i),drop=FALSE],
                                    2, sum )
            }
        }
    } else{
        ic <- list( positive=array( NA, c( maxorder+1, L ) ),
                    negative=array( NA, c( maxorder+1, L ) ) )
        for( i in 1:( maxorder+1 ) ){
            if( !( is.null( kmerContributions[[i]] ) ) ){
                matr <- kmerContributions[[i]][,-(1:i),drop=FALSE]
                matr[ matr < 0 ] <- 0
                ic[["positive"]][i,i:L] <- apply( matr, 2, sum )

                matr <- kmerContributions[[i]][,-(1:i),drop=FALSE]
                matr[ matr > 0 ] <- 0
                ic[["negative"]][i,i:L] <- apply( matr, 2, sum )
            }
        }
    }    
    return( ic )
}

# calculate the position-specific contribution of BaMM k-mers to the information
# content from the model probabilities <probs> and conditional probabilities
# <conds> for the orders <order> using background frequencies <freqs> for the
# 0'th order. The information content is computed by using the logarithm to the
# base <base>. Choose between RNA (<rna>: TRUE) and DNA (<rna>: FALSE) alphabet.

kmerContributionsToInformationContent <- function( probs, conds, order,
                                                   freqs=rep( .25 ,4 ), base=2,
                                                   rna=FALSE ){
    MAXORDER <- 5

    # model formats
    L <- unique( c( sapply( probs, nrow ), sapply( conds, nrow ) ) )
    maxorder <- unique( length( probs ), length( conds ) )-1
    if( length( L ) > 1 || length( maxorder ) > 1 ){
        stop( "probs and conds formats differ" )
    }
    
    # maximum model order
    if( any( order > maxorder ) ){
        missing.order <- order[order > maxorder]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " missing", sep="" ) )
    }

    # maximum code order
    if( any( order > MAXORDER ) ){
        missing.order <- order[order > MAXORDER]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " not implemented yet", sep="" ) )
    }

    # determine maximum order
    maxorder <- min( maxorder, MAXORDER )

    kmerContributions.list <- list()
    if( rna ){
        nucleotides <- c( "A", "C", "G", "U" )
    } else{
        nucleotides <- c( "A", "C", "G", "T" )
    }

    if( maxorder >= 0 && 0 %in% order ){

        kmerContributions <- list()
        for( a in 1:4 ){
            kmerContributions <- c( kmerContributions, list(
                                    probs[[1]][1:L,a] * log(
                                    probs[[1]][1:L,a] /
                                    freqs[a], base ) ) )
        }

        A <- rep( nucleotides , 1 )

        kmerContributions <- cbind( "i"=A,
                                    data.frame( matrix( unlist(
                                    kmerContributions ), nrow=4, byrow=TRUE )
                                    ), stringsAsFactors=FALSE )

        colnames( kmerContributions )[ -1 ] <- 1:L
        rownames( kmerContributions ) <- 1:4

        kmerContributions.list[[1]] <- kmerContributions
    }

    if( maxorder >= 1 && 1 %in% order ){

    kmerContributions <- list()

        for( b in 1:4 ){
            for( a in 1:4 ){
                kmerContributions <- c( kmerContributions, list(
                                        probs[[2]][2:L,a,b] * log(
                                        conds[[2]][2:L,a,b] /
                                        conds[[1]][2:L,b], base ) ) )
            }
        }

        A <- rep( nucleotides, 4 )
        B <- rep( as.character( matrix( nucleotides, nrow=4, ncol=4, byrow=TRUE
                  ) ), 1 )

        kmerContributions <- cbind( "i-1"=A, "i"=B,
                                    data.frame( matrix( unlist(
                                    kmerContributions ), nrow=16, byrow=TRUE )
                                    ), stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2] ), ]

        colnames( kmerContributions )[ -( 1:2 ) ] <- 2:L
        rownames( kmerContributions ) <- 1:16

        kmerContributions.list[[2]] <- kmerContributions
    }

    if( maxorder >= 2 && 2 %in% order ){
    
        kmerContributions <- list()
        for( c in 1:4 ){
            for( b in 1:4 ){
                for( a in 1:4 ){
                    kmerContributions <- c( kmerContributions, list(
                                            probs[[3]][3:L,a,b,c] * log(
                                            conds[[3]][3:L,a,b,c] /
                                            conds[[2]][3:L,b,c], base ) ) )
                }
            }
        }

        A <- rep( nucleotides, 16 )
        B <- rep( as.character( matrix( nucleotides, nrow= 4, ncol=4, byrow=TRUE
                  ) ), 4 )
        C <- rep( as.character( matrix( nucleotides, nrow=16, ncol=4, byrow=TRUE
                  ) ), 1 )

        kmerContributions <- cbind( "i-2"=A, "i-1"=B, "i"=C,
                                    data.frame( matrix( unlist(
                                    kmerContributions ), nrow=64, byrow=TRUE )
                                    ), stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3] ), ]

        colnames( kmerContributions )[ -( 1:3 ) ] <- 3:L
        rownames( kmerContributions ) <- 1:64

        kmerContributions.list[[3]] <- kmerContributions
    }

    if( maxorder >= 3 && 3 %in% order ){
    
        kmerContributions <- list()

        for( d in 1:4 ){
            for( c in 1:4 ){
                for( b in 1:4 ){
                    for( a in 1:4 ){
                        kmerContributions <- c( kmerContributions, list(
                                                probs[[4]][4:L,a,b,c,d] * log(
                                                conds[[4]][4:L,a,b,c,d] /
                                                conds[[3]][4:L,b,c,d], base ) ) )
                    }
                }
            }
        }

        A <- rep( nucleotides, 64 )
        B <- rep( as.character( matrix( nucleotides, nrow= 4, ncol=4, byrow=TRUE
                  ) ), 16 )
        C <- rep( as.character( matrix( nucleotides, nrow=16, ncol=4, byrow=TRUE
                  ) ),  4 )
        D <- rep( as.character( matrix( nucleotides, nrow=64, ncol=4, byrow=TRUE
                  ) ),  1 )

        kmerContributions <- cbind( "i-3"=A, "i-2"=B, "i-1"=C, "i"=D,
                                    data.frame( matrix( unlist(
                                    kmerContributions ), nrow=256, byrow=TRUE )
                                    ), stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3],
                                                kmerContributions[,4] ), ]

        colnames( kmerContributions )[ -( 1:4 ) ] <- 4:L
        rownames( kmerContributions ) <- 1:256

        kmerContributions.list[[4]] <- kmerContributions
    }

    if( maxorder >= 4 && 4 %in% order ){

        kmerContributions <- list()

        for( e in 1:4 ){
            for( d in 1:4 ){
                for( c in 1:4 ){
                    for( b in 1:4 ){
                        for( a in 1:4 ){
                            kmerContributions <- c( kmerContributions, list(
                                                    probs[[5]][5:L,a,b,c,d,e] * log(
                                                    conds[[5]][5:L,a,b,c,d,e] /
                                                    conds[[4]][5:L,b,c,d,e], base ) ) )
                        }
                    }
                }
            }
        }

        A <- rep( nucleotides, 256 )
        B <- rep( as.character( matrix( nucleotides, nrow=  4, ncol=4, byrow=TRUE
                  ) ), 64 )
        C <- rep( as.character( matrix( nucleotides, nrow= 16, ncol=4, byrow=TRUE
                  ) ), 16 )
        D <- rep( as.character( matrix( nucleotides, nrow= 64, ncol=4, byrow=TRUE
                  ) ),  4 )
        E <- rep( as.character( matrix( nucleotides, nrow=256, ncol=4, byrow=TRUE
                  ) ),  1 )

        kmerContributions <- cbind( "i-4"=A, "i-3"=B, "i-2"=C, "i-1"=D, "i"=E,
                                    data.frame( matrix( unlist(
                                    kmerContributions ), nrow=1024, byrow=TRUE )
                                    ), stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3],
                                                kmerContributions[,4],
                                                kmerContributions[,5] ), ]

        colnames( kmerContributions )[ -( 1:5 ) ] <- 5:L
        rownames( kmerContributions ) <- 1:1024

        kmerContributions.list[[5]] <- kmerContributions
    }

    if( maxorder >= 5 && 5 %in% order ){

        kmerContributions <- list()

        for( f in 1:4 ){
            for( e in 1:4 ){
                for( d in 1:4 ){
                    for( c in 1:4 ){
                        for( b in 1:4 ){
                            for( a in 1:4 ){
                                kmerContributions <- c( kmerContributions, list(
                                                        probs[[6]][6:L,a,b,c,d,e,f] * log(
                                                        conds[[6]][6:L,a,b,c,d,e,f] /
                                                        conds[[5]][6:L,b,c,d,e,f], base ) ) )
                            }
                        }
                    }
                }
            }
        }

        A <- rep( nucleotides, 1024 )
        B <- rep( as.character( matrix( nucleotides, nrow=   4, ncol=4,
                  byrow=TRUE ) ), 256 )
        C <- rep( as.character( matrix( nucleotides, nrow=  16, ncol=4,
                  byrow=TRUE ) ),  64 )
        D <- rep( as.character( matrix( nucleotides, nrow=  64, ncol=4,
                  byrow=TRUE ) ),  16 )
        E <- rep( as.character( matrix( nucleotides, nrow= 256, ncol=4,
                  byrow=TRUE ) ),   4 )
        F <- rep( as.character( matrix( nucleotides, nrow=1024, ncol=4,
                  byrow=TRUE ) ),   1 )

        kmerContributions <- cbind( "i-5"=A, "i-4"=B, "i-3"=C, "i-2"=D, "i-1"=E, "i"=F,
                                    data.frame( matrix( unlist(
                                    kmerContributions ), nrow=4096, byrow=TRUE
                                    ) ), stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3],
                                                kmerContributions[,4],
                                                kmerContributions[,5],
                                                kmerContributions[,6] ), ]

        colnames( kmerContributions )[ -( 1:6 ) ] <- 6:L
        rownames( kmerContributions ) <- 1:4096

        kmerContributions.list[[6]] <- kmerContributions
    }

    return( kmerContributions.list )
}

# comments in wrapper
hoSeqLogo <- function( filename, order, base=2, rna=FALSE, useFreqs=FALSE,
                       icColumnScale=TRUE, icLetterScale=TRUE, alpha=.3,
                       plot.xaxis=TRUE, plot.xlab=TRUE, plot.border.labels=TRUE,
                       xanchor=NULL, xanchor.labels=TRUE, plot.xanchor=TRUE,
                       plot.yaxis=TRUE, plot.ylab=TRUE, ylim=NULL, yat=NULL,
                       yaxis.vjust=.5, cex=2.7, lwd=4, probs.ending="probs",
                       conds.ending="conds", freqs.ending="freqs" ){

    if( !icColumnScale && icLetterScale ){
        stop( "icLetterScale and not icColumnScale is not implemented" )
    }

    # read in or set default background frequencies
    if( useFreqs ){
        file.freqs <- paste( filename, freqs.ending, sep="." )
        freqs <- scan( file.freqs, quiet=TRUE )
    } else{
        freqs <- rep( .25, 4 )
    }

    # read in conditional probabilities
    file.conds <- paste( filename, conds.ending, sep="." )
    conds <- readBaMM( file.conds )

    # read in probabilities
    file.probs <- paste( filename, probs.ending, sep="." )
    probs <- readBaMM( file.probs )

    MAXORDER <- 5

    # model formats
    L <- unique( c( sapply( probs, nrow ), sapply( conds, nrow ) ) )
    maxorder <- unique( length( probs ), length( conds ) ) - 1
    if( length( L ) > 1 || length( maxorder ) > 1 ){
        stop( paste( file.probs, " and ", file.conds, " formats differ", sep=""
              ) )
    }
    
    # maximum model order
    if( any( order > maxorder ) ){
        missing.order <- order[order > maxorder]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " missing.", sep="" ) )
    }

    # maximum code order
    if( any( order > MAXORDER ) ){
        missing.order <- order[order > MAXORDER]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " not implemented yet", sep="" ) )
    }

    # determine maximum order
    maxorder <- min( maxorder, MAXORDER, max( order ) )

    if( icColumnScale || icLetterScale ){
        # calculate the (position-specific) contribution of BaMM k-mers to the
        # information content
        kmerContributions <- kmerContributionsToInformationContent( probs, conds, order, freqs, base, rna )
    }

    wts <- rep( 1, MAXORDER+1 )
    wt.gps <- c( 0, rep( .4, MAXORDER ) )
    ht.gps <- c( .01, .005, rep( .001, MAXORDER-1 ) )

    if( rna ){
        chars <- c( "A", "C", "G", "U" )
    } else{
        chars <- c( "A", "C", "G", "T" )
    }

    letters.list <- list()

    if( icLetterScale ){

        if( is.null( ylim ) ){
            ylim <- c( -2, 2 )
        }
        ylab <- "Information content"

        # calculate the (position-specific) information content of the BaMM from
        # k-mer contributions
#         facs <- informationContent.kmerContributions( kmerContributions,
#                                                       unsigned=FALSE )

        k <- 0
        while( maxorder >= k ){

            if( !( k %in% order ) ){
                k <- k + 1
                next # skip this order's sequence logo
            }

            k1 <- k+1
            columns <- kmerContributions[[k1]][,-(1:k1),drop=FALSE]

            if( any( is.nan( as.matrix( columns ) ) ) ||
                any( is.infinite( as.matrix( columns ) ) ) ){

                warning( paste( "NaN and/or -Inf ", k1, "-mer contributions to",
                " the information content. Plot ", k, ifelse( k %in% 1:3,
                switch( k, "st", "nd", "rd" ), "th" ), "-order ",
                "sequence logo by using probability scale", sep="" ) )

                k <- k + 1
                next # skip this order's sequence logo
            }

            # p(o)s(iti)v(e)
            cols.psv <- columns
            cols.psv[ cols.psv < 0 ] <- 0

            # n(e)g(ati)v(e)
            cols.ngv <- columns
            cols.ngv[ cols.ngv > 0 ] <- 0
            cols.ngv <- ( -1 ) * cols.ngv

            wt <- wts[k1]
            wt.gp <- wt.gps[k1]

            chars.matrix <- as.matrix( kmerContributions[[k1]][,1:k1,drop=FALSE]
                                       )

            x.pos <- 0
            letters <- list( x=NULL, y=NULL, id=NULL, fill=NULL, col=NULL )

            for( j in 1:(L-k) ){

                ht.gp <- ht.gps[k1] # reset letter gaps

                scale.factor <- 1-( length( which( cols.psv[,j] > 0 ) ) * ht.gp
                                    )
                # scale letters
                hts <- scale.factor * cols.psv[,j]
                # scale letter gaps
                ht.gp <- ht.gp * sum( cols.psv[,j] )
                letterOrder <- order( hts )

                y.pos <- ht.gp

                for( i in 1:nrow( chars.matrix ) ){
                    ht <- hts[ letterOrder[i] ]
                    if( ht > 0 ){
                        k1mer <- chars.matrix[ letterOrder[i], ]
                        for( l in 1:k1 ){
                            letters <- addLetter( letters, k1mer[l],
                                                  x.pos+( ( l-1 )*wt ), y.pos,
                                                  ht, wt )
                        }
                        y.pos <- y.pos + ht + ht.gp
                    }
                }

                ht.gp <- ht.gps[k1] # reset letter gaps

                scale.factor <- 1-( length( which( cols.ngv[,j] > 0 ) ) * ht.gp
                                    )
                # scale letters
                hts <- scale.factor * cols.ngv[,j]
                # scale letter gaps
                ht.gp <- ht.gp * sum( cols.ngv[,j] )
                letterOrder <- order( hts )

                y.pos <- -ht.gp
                for( i in 1:nrow( chars.matrix ) ){
                    ht <- hts[ letterOrder[i] ]
                    if( ht > 0 ){
                        y.pos <- y.pos - ht
                        k1mer <- chars.matrix[ letterOrder[i], ]
                        for( l in 1:k1 ){
                            letters <- addLetter( letters, k1mer[l],
                                                  x.pos+( ( l-1 )*wt ), y.pos,
                                                  ht, wt )
                        }
                        y.pos <- y.pos - ht.gp
                    }
                }
                x.pos <- x.pos + ( k1*wt ) + wt.gp
            }
            letters.list[[k1]] <- letters
            k <- k + 1
        }
    } else{

        if( icColumnScale ){

            if( is.null( ylim ) ){
                ylim <- c( 0, 2 )
            }
            ylab <- "Information content"

            # calculate the (position-specific) information content of the BaMM
            # from k-mer contributions
            facs <- informationContent.kmerContributions( kmerContributions,
                                                          unsigned=TRUE )
        } else{

            if( is.null( ylim ) ){
                ylim <- c( 0, 1 )
            }
            ylab <- "Probability"

            # (position-specific) k-mer probabilities add up to one
            facs <- matrix( 1, nrow=maxorder+1, L )
        }

        k <- 0
        while( maxorder >= k ){

            if( !( k %in% order ) ){
                k <- k + 1
                next # skip this order's sequence logo
            }

            k1 <- k+1
            columns <- probs[[k1]]

            if( icColumnScale ){
                if( any( is.nan( facs[k1,k1:L] ) ) ||
                    any( is.infinite( facs[k1,k1:L] ) ) ){

                    warning( paste( "NaN and/or -Inf ", k1, "-mer ",
                             "contributions to the information content. Plot ",
                             k, ifelse( k %in% 1:3, switch( k, "st", "nd", "rd"
                             ), "th" ), "-order sequence logo by using ",
                             "probability scale.", sep="" ) )

                    k <- k + 1
                    next # skip this order's sequence logo
                }
            }

            wt <- wts[k1]

            wt.gp <- wt.gps[k1]

            chars.matrix <- matrix( rep( chars, 4^k ), ncol=1 )
            if( k > 0 ){
                for( i in 1:k ){
                    chars.matrix <- cbind( chars.matrix, rep( as.character(
                                           matrix( chars, nrow=4^i, ncol=4,
                                           byrow=TRUE ) ), 4^(k-i) ) )
                }
            }

            x.pos <- 0
            letters <- list( x=NULL, y=NULL, id=NULL, fill=NULL, col=NULL )

            for( j in k1:L ){

                ht.gp <- ht.gps[k1] # reset letter gaps

                column <- switch( k1, columns[j,], columns[j,,], columns[j,,,],
                                  columns[j,,,,], columns[j,,,,,],
                                  columns[j,,,,,,] )

                scale.factor <- 1-( ( length( which( column > 0 ) )-1 ) * ht.gp
                                    )
                # scale letters
                hts <- scale.factor * column * facs[k1,j]
                # scale letter gaps
                ht.gp <- ht.gp * facs[k1,j]
                letterOrder <- order( hts )

                y.pos <- 0
                for( i in 1:nrow( chars.matrix ) ){
                    ht <- hts[ letterOrder[i] ]
                    if( ht > 0 ){
                        k1mer <- chars.matrix[ letterOrder[i], ]
                        for( l in 1:k1 ){
                            letters <- addLetter( letters, k1mer[l],
                                                  x.pos+( ( l-1 )*wt ), y.pos,
                                                  ht, wt )
                        }
                        y.pos <- y.pos + ht + ht.gp
                    }
                }
                x.pos <- x.pos + ( k1*wt ) + wt.gp
            }

            letters.list[[k1]] <- letters
            k <- k + 1
        }
    }

    for( i in 1:length( letters.list ) ){

        letters <- letters.list[[i]]

        if( unique( is.null( unlist( letters ) ) ) ){
            if( ( i-1 ) %in% order ){
                cat( "Check ", i, "-mer contributions to the information ",
                     "content\n", sep="" )
                print( kmerContributions[[i]] )
            }
            next # no letter coordinates
        }

        unit <- i * wts[i] # width of one i-mer
        units <- ( L-(i-1) ) * unit # width of all i-mers

        gap <- wt.gps[i] # width of one gap
        gaps <- ( L-i ) * gap # width of all gaps

        offset <- unit / 2

#         grid.newpage()

        yaxis.transparent <- FALSE
        if( is.character( plot.yaxis ) ){
            if( plot.yaxis == "transparent" ){
                plot.yaxis <- TRUE
                yaxis.transparent <- TRUE
            }
        }

        bottomMargin = ifelse( plot.xaxis, 8, 2 )
        leftMargin = ifelse( plot.yaxis, 8, 2 )
        
        topMargin = 6
#         topMargin = 12

        pushViewport( plotViewport( c( bottomMargin, leftMargin, topMargin, 2 )
                      ) )
        pushViewport( dataViewport( 0:(units+gaps), c( ylim[1], ylim[2] ),
                      name="vp1" ) )

#         grid.polygon( x=unit( c( 0, units+gaps+1, units+gaps+1, 0 ), "native" ),
#                       y=unit( c( 0, 0, 2, 2 ), "native" ), gp=gpar(
#                       fill="transparent", col="transparent" ) )
        grid.polygon( x=unit( c( 0, units+gaps+1, units+gaps+1, 0 ), "native" ),
                      y=unit( c( 0, 0, ylim[2], ylim[2] ), "native" ), gp=gpar(
                      fill="transparent", col="transparent" ) )
        grid.polygon( x=unit( letters$x,"native" ), y=unit( letters$y, "native"
                      ), id=letters$id, gp=gpar( fill=letters$fill,
                      col=letters$col, lwd=1 ) )

        if( icLetterScale ){
            # superimpose white transparent polygon
            if( unit == units ){
                xunit <- unit( c( -.2, units+gaps, units+gaps, -.2 ),
                               "native" )
            } else{
                xunit <- unit( c( -.2, units+gaps+1, units+gaps+1, -.2 ),
                               "native" )
            }
            grid.polygon( x=xunit, y=unit( c( 0, 0, ylim[1], ylim[1] ),
                          "native" ), gp=gpar( fill="#FFFFFF", col="#FFFFFF",
                          alpha=1-alpha ) )
        }

        if( plot.xaxis ){
            if( unit == units ){
                if( is.null( xanchor ) ){
                    grid.xaxis( at=offset, label=1, edits=gEdit( gPath="labels",
                                vjust=0 ), gp=gpar( cex=cex, lwd=lwd ) )
                } else{
                    grid.xaxis( at=offset, label=xanchor.labels, edits=gEdit(
                                gPath="labels", vjust=0 ), gp=gpar( cex=cex,
                                lwd=lwd ) )
                }
            } else{
                if( is.null( xanchor ) ){
                    if( plot.border.labels ){
                        grid.xaxis( at=c( offset, units+gaps-offset ), label=c(
                                    i, L ), gp=gpar( cex=cex, lwd=lwd ) )
                    } else{
                        grid.xaxis( at=seq( offset, units+gaps-offset, unit+gap
                                    ), label=i:L, gp=gpar( cex=cex, lwd=lwd ) )
                    }
                } else{
                    if( is.numeric( xanchor.labels ) ){
                        if( plot.xanchor ){
                            if( length( xanchor.labels ) == 1 ){
                                grid.xaxis( at=c( offset, ( ( xanchor-1-( i-1 )
                                            )*unit )+( ( xanchor-1-( i-1 ) )*gap
                                            )+offset, units+gaps-offset ),
                                            label=c( xanchor.labels-( (
                                            xanchor-1 )-( i-1 ) ),
                                            xanchor.labels, xanchor.labels+(
                                            xanchor-1 ) ), edits=gEdit(
                                            gPath="labels", vjust=0 ), gp=gpar(
                                            cex=cex, lwd=lwd ) )
                            } else if ( length( xanchor.labels == 3 ) ){
                                grid.xaxis( at=c( offset, ( ( xanchor-1-( i-1 )
                                            )*unit )+( ( xanchor-1-( i-1 ) )*gap
                                            )+offset, units+gaps-offset ),
                                            label=xanchor.labels, edits=gEdit(
                                            gPath="labels", vjust=0 ), gp=gpar(
                                            cex=cex, lwd=lwd ) )
                            } else{
                                stop( "argument xanchor.labels must be a vector of length 1 or 3" )
                            }
                        } else{
#                             grid.xaxis( at=c( offset, units+gaps-offset ),
#                                         label=c( xanchor.labels-( ( xanchor-1
#                                         )-( i-1 ) ), xanchor.labels+( xanchor-1
#                                         ) ), edits=gEdit( gPath="labels",
#                                         vjust=0 ), gp=gpar( cex=cex, lwd=lwd ) )
                            grid.xaxis( at=c( offset, units+gaps-offset ),
                                        label=c( xanchor.labels-( ( xanchor-1
                                        )-( i-1-order ) ), xanchor.labels+(
                                        xanchor-1 ) ), edits=gEdit(
                                        gPath="labels", vjust=0 ), gp=gpar(
                                        cex=cex, lwd=lwd ) )
                        }
                    } else{
                        if( length( xanchor.labels ) == 1 ){
                            grid.xaxis( at=c( offset, ( ( xanchor-1-( i-1 )
                                        )*unit )+( ( xanchor-1-( i-1 ) )*gap
                                        )+offset, units+gaps-offset ), label=c(
                                        -( xanchor-1-( i-1 ) ), xanchor.labels,
                                        units/unit-( xanchor-( i-1 ) ) ),
                                        edits=gEdit( gPath="labels", vjust=0 ),
                                        gp=gpar( cex=cex, lwd=lwd ) )
                        } else if ( length( xanchor.labels ) == 3 ){
                            grid.xaxis( at=c( offset, ( ( xanchor-1-( i-1 )
                                        )*unit )+( ( xanchor-1-( i-1 ) )*gap
                                        )+offset, units+gaps-offset ),
                                        label=xanchor.labels, edits=gEdit(
                                        gPath="labels", vjust=0 ), gp=gpar(
                                        cex=cex, lwd=lwd ) )
                        } else{
                            stop( "argument xanchor.labels must be a vector of length 1 or 3" )
                        }
                    }
                }
            }

            if( plot.xlab ){
                grid.text( "Model position [nt]", y=unit( -2.5, "lines" ),
                           gp=gpar( cex=cex ) )
            }
        }

        if( plot.yaxis ){
            if( !( yaxis.transparent ) ){
                if( is.null( yat ) ){
                    grid.yaxis( gp=gpar( cex=cex, lwd=lwd ) )
                    if( plot.ylab ){
                        grid.text( ylab, x=unit( -2.5, "lines" ), rot=90,
                                   gp=gpar( cex=cex ) )
                    }
                } else{
                    grid.yaxis( at=yat, label=TRUE, edits=gEdit( gPath="labels",
                                hjust=.5, vjust=0, rot=90 ), gp=gpar( cex=cex,
                                lwd=lwd ) )
                    if( plot.ylab ){
                        grid.text( ylab, x=unit( -2.5, "lines" ),
                                   hjust=yaxis.vjust, rot=90, gp=gpar( cex=cex )
                                   )
                    }
                }
            }
        }
        popViewport()
        popViewport()

        par( ask=FALSE )
    }
}

# calculate the higher-order information content of the homogeneous BaMM from
# k-mer contributions calculated by <bgKmerContributionsToInformationContent>.
# Discriminate between positive and negative contributions (<unsigned>: FALSE)
# or not (<unsigned: TRUE).
informationContent.bgKmerContributions <- function( bgKmerContributions,
                                                    unsigned=TRUE ){

    maxorder <- length( bgKmerContributions )-1
    # n(ot) n(ull) indices
    nnindices <- which( !( sapply( bgKmerContributions, is.null ) ) )
    L <- unique( sapply( bgKmerContributions[ nnindices ], ncol ) - ( 1:(
                 maxorder+1 ) ) )

    if( length( L ) > 1 || L != 1 ){
        stop( "bgKmerContributions format error" )
    }

    if( unsigned ){
        ic <- rep( NA, maxorder+1 )

        for( i in 1:( maxorder+1 ) ){
            if( !( is.null( bgKmerContributions[[i]] ) ) ){
                ic[i] <- sum( bgKmerContributions[[i]][,-(1:i)] )
            }
        }
    } else{
        ic <- list( positive=rep( NA, maxorder+1 ),
                    negative=rep( NA, maxorder+1 ) )

        for( i in 1:( maxorder+1 ) ){
            if( !( is.null( bgKmerContributions[[i]] ) ) ){

                vec <- bgKmerContributions[[i]][,-(1:i)]
                vec[ vec < 0 ] <- 0
                ic[["positive"]][i] <- sum( vec )

                vec <- bgKmerContributions[[i]][,-(1:i)]
                vec[ vec > 0 ] <- 0
                ic[["negative"]][i] <- sum( vec )
            }
        }
    }
    return( ic )
}

# calculate the contribution of homogeneous BaMM k-mers to the information
# content from the model probabilities <probs> and conditional probabilities
# <conds> for the orders <order> using background frequencies <freqs> for the
# 0'th order. The information content is computed by using the logarithm to the
# base <base>.
bgKmerContributionsToInformationContent <- function( probs, conds, order,
                                                     freqs=rep( .25 ,4 ),
                                                     base=2, rna=FALSE ){
    MAXORDER <- 6

    # model formats
    L <- unique( c( sapply( probs, nrow ), sapply( conds, nrow ) ) )
    maxorder <- unique( length( probs ), length( conds ) )-1
    if( length( L ) > 1 || length( maxorder ) > 1 ){
        stop( "probs and conds formats differ" )
    }
    if( L != 1 ){
        stop( "BaMM format error" )
    }
    
    # maximum model order
    if( any( order > maxorder ) ){
        missing.order <- order[order > maxorder]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " missing.", sep="" ) )
    }

    # maximum code order
    if( any( order > MAXORDER ) ){
        missing.order <- order[order > MAXORDER]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " not implemented yet.", sep="" ) )
    }

    # determine maximum order
    maxorder <- min( maxorder, MAXORDER )

    kmerContributions.list <- list()
    if( rna ){
        nucleotides <- c( "A", "C", "G", "U" )
    } else{
        nucleotides <- c( "A", "C", "G", "T" )
    }

    if( maxorder >= 0 && 0 %in% order ){
    
        kmerContributions <- double( 4 )

        k <- 0
        for( a in 1:4 ){
            k <- k + 1
            kmerContributions[k] <- probs[[1]][,a] * log(
                                    probs[[1]][,a] /
                                    freqs[a], base )
        }

        A <- rep( nucleotides , 1 )

        kmerContributions <- data.frame( "i"=A,
                                         kmerContributions,
                                         stringsAsFactors=FALSE )

        colnames( kmerContributions )[2] <- 1
        rownames( kmerContributions ) <- 1:4

        kmerContributions.list[[1]] <- kmerContributions
    }

    if( maxorder >= 1 && 1 %in% order ){
    
        kmerContributions <- double( 16 )

        k <- 0
        for( b in 1:4 ){
            for( a in 1:4 ){
                k <- k + 1
                kmerContributions[k] <- probs[[2]][,a,b] * log(
                                        conds[[2]][,a,b] /
                                        conds[[1]][,b], base )
            }
        }

        A <- rep( nucleotides, 4 )
        B <- rep( as.character( matrix( nucleotides, nrow=4, ncol=4, byrow=TRUE
                  ) ), 1 )

        kmerContributions <- data.frame( "i-1"=A, "i"=B,
                                         kmerContributions,
                                         stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2] ), ]

        colnames( kmerContributions )[3] <- 1
        rownames( kmerContributions ) <- 1:16

        kmerContributions.list[[2]] <- kmerContributions
    }

    if( maxorder >= 2 && 2 %in% order ){
    
        kmerContributions <- double( 64 )

        k <- 0
        for( c in 1:4 ){
            for( b in 1:4 ){
                for( a in 1:4 ){
                    k <- k + 1
                    kmerContributions[k] <- probs[[3]][,a,b,c] * log(
                                            conds[[3]][,a,b,c] /
                                            conds[[2]][,b,c], base )
                }
            }
        }

        A <- rep( nucleotides, 16 )
        B <- rep( as.character( matrix( nucleotides, nrow= 4, ncol=4, byrow=TRUE
                  ) ), 4 )
        C <- rep( as.character( matrix( nucleotides, nrow=16, ncol=4, byrow=TRUE
                  ) ), 1 )

        kmerContributions <- data.frame( "i-2"=A, "i-1"=B, "i"=C,
                                         kmerContributions,
                                         stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3] ), ]

        colnames( kmerContributions )[4] <- 1
        rownames( kmerContributions ) <- 1:64

        kmerContributions.list[[3]] <- kmerContributions
    }

    if( maxorder >= 3 && 3 %in% order ){
    
        kmerContributions <- double( 256 )

        k <- 0
        for( d in 1:4 ){
            for( c in 1:4 ){
                for( b in 1:4 ){
                    for( a in 1:4 ){
                        k <- k + 1
                        kmerContributions[k] <- probs[[4]][,a,b,c,d] * log(
                                                conds[[4]][,a,b,c,d] /
                                                conds[[3]][,b,c,d], base )
                    }
                }
            }
        }

        A <- rep( nucleotides, 64 )
        B <- rep( as.character( matrix( nucleotides, nrow= 4, ncol=4, byrow=TRUE
                  ) ), 16 )
        C <- rep( as.character( matrix( nucleotides, nrow=16, ncol=4, byrow=TRUE
                  ) ),  4 )
        D <- rep( as.character( matrix( nucleotides, nrow=64, ncol=4, byrow=TRUE
                  ) ),  1 )

        kmerContributions <- data.frame( "i-3"=A, "i-2"=B, "i-1"=C, "i"=D,
                                         kmerContributions,
                                         stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3],
                                                kmerContributions[,4] ), ]

        colnames( kmerContributions )[5] <- 1
        rownames( kmerContributions ) <- 1:256

        kmerContributions.list[[4]] <- kmerContributions
    }
    
    if( maxorder >= 4 && 4 %in% order ){

        kmerContributions <- double( 1024 )

        k <- 0
        for( e in 1:4 ){
            for( d in 1:4 ){
                for( c in 1:4 ){
                    for( b in 1:4 ){
                        for( a in 1:4 ){
                            k <- k + 1
                            kmerContributions[k] <- probs[[5]][,a,b,c,d,e] * log(
                                                    conds[[5]][,a,b,c,d,e] /
                                                    conds[[4]][,b,c,d,e], base )
                        }
                    }
                }
            }
        }

        A <- rep( nucleotides,   256 )
        B <- rep( as.character( matrix( nucleotides, nrow=  4, ncol=4,
                  byrow=TRUE ) ), 64 )
        C <- rep( as.character( matrix( nucleotides, nrow= 16, ncol=4,
                  byrow=TRUE ) ), 16 )
        D <- rep( as.character( matrix( nucleotides, nrow= 64, ncol=4,
                  byrow=TRUE ) ),  4 )
        E <- rep( as.character( matrix( nucleotides, nrow=256, ncol=4,
                  byrow=TRUE ) ),  1 )

        kmerContributions <- data.frame( "i-4"=A, "i-3"=B, "i-2"=C, "i-1"=D, "i"=E,
                                         kmerContributions,
                                         stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3],
                                                kmerContributions[,4],
                                                kmerContributions[,5] ), ]

        colnames( kmerContributions )[6] <- 1
        rownames( kmerContributions ) <- 1:1024

        kmerContributions.list[[5]] <- kmerContributions
    }

    if( maxorder >= 5 && 5 %in% order ){

        kmerContributions <- double( 4096 )

        k <- 0
        for( f in 1:4 ){
            for( e in 1:4 ){
                for( d in 1:4 ){
                    for( c in 1:4 ){
                        for( b in 1:4 ){
                            for( a in 1:4 ){
                                k <- k + 1
                                kmerContributions[k] <- probs[[6]][,a,b,c,d,e,f] * log(
                                                        conds[[6]][,a,b,c,d,e,f] /
                                                        conds[[5]][,b,c,d,e,f], base )
                            }
                        }
                    }
                }
            }
        }

        A <- rep( nucleotides,   1024 )
        B <- rep( as.character( matrix( nucleotides, nrow=   4, ncol=4,
                  byrow=TRUE ) ), 256 )
        C <- rep( as.character( matrix( nucleotides, nrow=  16, ncol=4,
                  byrow=TRUE ) ),  64 )
        D <- rep( as.character( matrix( nucleotides, nrow=  64, ncol=4,
                  byrow=TRUE ) ),  16 )
        E <- rep( as.character( matrix( nucleotides, nrow= 256, ncol=4,
                  byrow=TRUE ) ),   4 )
        F <- rep( as.character( matrix( nucleotides, nrow=1024, ncol=4,
                  byrow=TRUE ) ),   1 )

        kmerContributions <- data.frame( "i-5"=A, "i-4"=B, "i-3"=C, "i-2"=D, "i-1"=E, "i"=F,
                                         kmerContributions,
                                         stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3],
                                                kmerContributions[,4],
                                                kmerContributions[,5],
                                                kmerContributions[,6] ), ]

        colnames( kmerContributions )[7] <- 1
        rownames( kmerContributions ) <- 1:4096

        kmerContributions.list[[6]] <- kmerContributions
    }

    if( maxorder >= 6 && 6 %in% order ){

        kmerContributions <- double( 16384 )

        k <- 0
        for( g in 1:4 ){
            for( f in 1:4 ){
                for( e in 1:4 ){
                    for( d in 1:4 ){
                        for( c in 1:4 ){
                            for( b in 1:4 ){
                                for( a in 1:4 ){
                                    k <- k + 1
                                    kmerContributions[k] <- probs[[7]][,a,b,c,d,e,f,g] * log(
                                                            conds[[7]][,a,b,c,d,e,f,g] /
                                                            conds[[6]][,b,c,d,e,f,g], base )
                                }
                            }
                        }
                    }
                }
            }
        }

        A <- rep( nucleotides,    4096 )
        B <- rep( as.character( matrix( nucleotides, nrow=   4, ncol=4,
                  byrow=TRUE ) ), 1024 )
        C <- rep( as.character( matrix( nucleotides, nrow=  16, ncol=4,
                  byrow=TRUE ) ),  256 )
        D <- rep( as.character( matrix( nucleotides, nrow=  64, ncol=4,
                  byrow=TRUE ) ),   64 )
        E <- rep( as.character( matrix( nucleotides, nrow= 256, ncol=4,
                  byrow=TRUE ) ),   16 )
        F <- rep( as.character( matrix( nucleotides, nrow=1024, ncol=4,
                  byrow=TRUE ) ),    4 )
        G <- rep( as.character( matrix( nucleotides, nrow=4096, ncol=4,
                  byrow=TRUE ) ),    1 )

        kmerContributions <- data.frame( "i-6"=A, "i-5"=B, "i-4"=C, "i-3"=D, "i-2"=E, "i-1"=F, "i"=G,
                                         kmerContributions,
                                         stringsAsFactors=FALSE )
        kmerContributions <- kmerContributions[ order(
                                                kmerContributions[,1],
                                                kmerContributions[,2],
                                                kmerContributions[,3],
                                                kmerContributions[,4],
                                                kmerContributions[,5],
                                                kmerContributions[,6],
                                                kmerContributions[,7] ), ]

        colnames( kmerContributions )[8] <- 1
        rownames( kmerContributions ) <- 1:16384

        kmerContributions.list[[7]] <- kmerContributions
    }

    return( kmerContributions.list )
}

# comments in wrapper
hoBgSeqLogo <- function( filename, order, base=2, rna=FALSE, icColumnScale=TRUE,
                         icLetterScale=TRUE, alpha=.3, plot.xaxis=TRUE,
                         plot.xlab=TRUE, plot.yaxis=TRUE, plot.ylab=TRUE,
                         ylim=NULL, yat=NULL, yaxis.vjust=.5, cex=2.7, lwd=4,
                         probs.ending="hbp", conds.ending="hbcp" ){

    if( !icColumnScale && icLetterScale ){
        stop( "icLetterScale and not icColumnScale is not implemented" )
    }

    # set default background frequencies
    freqs <- rep( .25, 4 )

    # read in conditional probabilities
    file.conds <- paste( filename, conds.ending, sep="." )
    conds <- readBaMM( file.conds )

    # read in probabilities
    file.probs <- paste( filename, probs.ending, sep="." )
    probs <- readBaMM( file.probs )

    MAXORDER <- 6

    # model formats
    L <- unique( c( sapply( probs, nrow ), sapply( conds, nrow ) ) )
    maxorder <- unique( length( probs ), length( conds ) ) - 1
    if( length( L ) > 1 || length( maxorder ) > 1 ){
        stop( paste( file.probs, " and ", file.conds, " formats differ", sep=""
              ) )
    }
    if( L != 1 ){
        stop( "BaMM format error" )
    }
   
    # maximum model order
    if( any( order > maxorder ) ){
        missing.order <- order[order > maxorder]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " missing", sep="" ) )
    }

    # maximum code order
    if( any( order > MAXORDER ) ){
        missing.order <- order[order > MAXORDER]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " not implemented yet", sep="" ) )
    }

    # determine maximum order
    maxorder <- min( maxorder, MAXORDER, max( order ) )

    if( icColumnScale || icLetterScale ){
        # calculate the contribution of homogeneous BaMM k-mers to the
        # information content
        kmerContributions <- bgKmerContributionsToInformationContent( probs, conds, order, freqs, base, rna )
    }

    wts <- rep( 1, MAXORDER+1 )
    ht.gps <- c( .01, .005, rep( .001, MAXORDER-1 ) )

    if( rna ){
        chars <- c( "A", "C", "G", "U" )
    } else{
        chars <- c( "A", "C", "G", "T" )
    }
    letters.list <- list()

    if( icLetterScale ){

        if( is.null( ylim ) ){
            ylim <- c( -2, 2 )
        }
        ylab <- "Information content"

        # calculate the information content of the homogeneous BaMM from k-mer
        # contributions
#         facs <- informationContent.bgKmerContributions( kmerContributions,
#                                                         unsigned=FALSE )

        k <- 0
        while( maxorder >= k ){

            if( !( k %in% order ) ){
                k <- k + 1
                next # skip this order's sequence logo
            }

            k1 <- k+1
            columns <- kmerContributions[[k1]][,-(1:k1)]

            if( any( is.nan( as.vector( columns ) ) ) ||
                any( is.infinite( as.vector( columns ) ) ) ){

                warning( paste( "NaN and/or -Inf ", k1, "-mer contributions to",
                         " the information content. Plot ", k, ifelse( k %in%
                         1:3, switch( k, "st", "nd", "rd" ), "th" ), "-order ",
                         "sequence logo by using probability scale.", sep="" ) )
                k <- k + 1
                next # skip this order's sequence logo
            }

            # p(o)s(iti)v(e)
            cols.psv <- columns
            cols.psv[ cols.psv < 0 ] <- 0

            # n(e)g(ati)v(e)
            cols.ngv <- columns
            cols.ngv[ cols.ngv > 0 ] <- 0
            cols.ngv <- ( -1 ) * cols.ngv

            wt <- wts[k1]
            ht.gp <- ht.gps[k1]

            chars.matrix <- as.matrix( kmerContributions[[k1]][,1:k1,drop=FALSE]
                                       )

            x.pos <- 0
            letters <- list( x=NULL, y=NULL, id=NULL, fill=NULL, col=NULL )

            scale.factor <- 1-( length( which( cols.psv > 0 ) ) * ht.gp )
            # scale letters
            hts <- scale.factor * cols.psv
            # scale letter gaps
            ht.gp <- ht.gp * sum( cols.psv )
            letterOrder <- order( hts )

            y.pos <- ht.gp
            for( i in 1:nrow( chars.matrix ) ){
                ht <- hts[ letterOrder[i] ]
                if( ht > 0 ){
                    k1mer <- chars.matrix[ letterOrder[i], ]
                    for( l in 1:k1 ){
                        letters <- addLetter( letters, k1mer[l],
                                              x.pos+( ( l-1 )*wt ), y.pos, ht,
                                              wt )
                    }
                    y.pos <- y.pos + ht + ht.gp
                }
            }

            ht.gp <- ht.gps[k1] # reset letter gaps

            scale.factor <- 1-( length( which( cols.ngv > 0 ) ) * ht.gp )
            # scale letters
            hts <- scale.factor * cols.ngv
            # scale letter gaps
            ht.gp <- ht.gp * sum( cols.ngv )
            letterOrder <- order( hts )

            y.pos <- -ht.gp
            for( i in 1:nrow( chars.matrix ) ){
                ht <- hts[ letterOrder[i] ]
                if( ht > 0 ){
                    y.pos <- y.pos - ht
                    k1mer <- chars.matrix[ letterOrder[i], ]
                    for( l in 1:k1 ){
                        letters <- addLetter( letters, k1mer[l],
                                              x.pos+( ( l-1 )*wt ), y.pos, ht,
                                              wt )
                    }
                    y.pos <- y.pos - ht.gp
                }
            }

            letters.list[[k1]] <- letters
            k <- k + 1
        }
    } else{

        if( icColumnScale ){

            if( is.null( ylim ) ){
                ylim <- c( 0, 2 )
            }
            ylab <- "Information content"

            # calculate the information content of the homogeneous BaMM from
            # k-mer contributions
            facs <- informationContent.bgKmerContributions( kmerContributions,
                                                            unsigned=TRUE )
        } else{

            if( is.null( ylim ) ){
                ylim <- c( 0, 1 )
            }
            ylab <- "Probability"

            # k-mer probabilities add up to one
            facs <- rep( 1, maxorder+1 )
        }

        k <- 0
        while( maxorder >= k ){

            if( !( k %in% order ) ){
                k <- k + 1
                next # skip this order's sequence logo
            }

            k1 <- k+1
            columns <- probs[[k1]]

            if( icColumnScale ){
                if( any( is.nan( facs[k1] ) ) ||
                    any( is.infinite( facs[k1] ) ) ){

                    warning( paste( "NaN and/or -Inf ", k1, "-mer ",
                             "contributions to the information content. Plot ",
                             k, ifelse( k %in% 1:3, switch( k, "st", "nd", "rd"
                             ), "th" ), "-order sequence logo by using ",
                             "probability scale.", sep="" ) )
                    k <- k + 1
                    next # skip this order's sequence logo
                }
            }

            wt <- wts[k1]
            ht.gp <- ht.gps[k1]

            chars.matrix <- matrix( rep( chars, 4^k ), ncol=1 )
            if( k > 0 ){
                for( i in 1:k ){
                    chars.matrix <- cbind( chars.matrix, rep( as.character(
                                           matrix( chars, nrow=4^i, ncol=4,
                                           byrow=TRUE ) ), 4^(k-i) ) )
                }
            }

            x.pos <- 0
            letters <- list( x=NULL, y=NULL, id=NULL, fill=NULL, col=NULL )

            column <- as.vector( columns )
            scale.factor <- 1-( ( length( which( column > 0 ) )-1 ) * ht.gp )
            # scale letters
            hts <- scale.factor * column * facs[k1]
            # scale letter gaps
            ht.gp <- ht.gp * facs[k1]
            letterOrder <- order( hts )

            y.pos <- 0
            for( i in 1:nrow( chars.matrix ) ){
                ht <- hts[ letterOrder[i] ]
                if( ht > 0 ){
                    k1mer <- chars.matrix[ letterOrder[i], ]
                    for( l in 1:k1 ){
                        letters <- addLetter( letters, k1mer[l],
                                              x.pos+( ( l-1 )*wt ), y.pos, ht,
                                              wt )
                    }
                    y.pos <- y.pos + ht + ht.gp
                }
            }

            letters.list[[k1]] <- letters
            k <- k + 1
        }
    }

    for( i in 1:length( letters.list ) ){

        letters <- letters.list[[i]]

        if( unique( is.null( unlist( letters ) ) ) ){
            if( ( i-1 ) %in% order ){
                cat( "Check ", i, "-mer contributions to the information ",
                "content\n", sep="" )
                print( kmerContributions[[i]] )
            }
            next # no letter coordinates
        }

        unit <- i * wts[i] # width of one i-mer

        grid.newpage()

        bottomMargin = ifelse( plot.xaxis, 8, 2 )
        leftMargin = ifelse( plot.yaxis, 8, 2 )

        pushViewport( plotViewport( c( bottomMargin, leftMargin, 2, 2 ) ) )
        pushViewport( dataViewport( 0:unit, c( ylim[1], ylim[2] ), name="vp1" )
                      )

        grid.polygon( x=unit( c( 0, unit+1, unit+1, 0 ), "native" ), y=unit( c(
                      0, 0, ylim[2], ylim[2] ), "native" ), gp=gpar(
                      fill="transparent", col="transparent" ) )
        grid.polygon( x=unit( letters$x,"native" ), y=unit( letters$y, "native"
                      ), id=letters$id, gp=gpar( fill=letters$fill,
                      col=letters$col, lwd=1 ) )

        if( icLetterScale ){
            # superimpose white transparent polygon
            grid.polygon( x=unit( c( -.2, unit+1, unit+1, -.2 ), "native" ),
                          y=unit( c( 0, 0, ylim[1], ylim[1] ), "native" ),
                          gp=gpar( fill="#FFFFFF", col="#FFFFFF", alpha=1-alpha
                          ) )
        }

        if( plot.xaxis ){
            grid.text( i-1, y=unit( -.3, "lines" ), gp=gpar( cex=cex, lwd=lwd )
                       )
            if( plot.xlab ){
                grid.text( "Order", y=unit( -1.5, "lines" ), gp=gpar(
                           cex=cex, lwd=lwd ) )
            }
        }
        if( plot.yaxis ){
            if( is.null( yat ) ){
                grid.yaxis( gp=gpar( cex=cex, lwd=lwd ) )
                if( plot.ylab ){
                    grid.text( ylab, x=unit( -2.5, "lines" ), rot=90, gp=gpar(
                               cex=cex, lwd=lwd ) )
                }
            } else{
                grid.yaxis( at=yat, label=TRUE, edits=gEdit( gPath="labels",
                            hjust=.5, vjust=0, rot=90 ), gp=gpar( cex=cex,
                            lwd=lwd ) )
                if( plot.ylab ){
                    grid.text( ylab, x=unit( -2.5, "lines" ), hjust=yaxis.vjust,
                               rot=90, gp=gpar( cex=cex ) )
                }
           }
        }
        popViewport()
        popViewport()

        par( ask=FALSE )
    }
}

hoContributionsToInformationContent <- function( filename, base=2, rna=FALSE,
                                                 useFreqs=FALSE,
                                                 probs.ending="probs",
                                                 conds.ending="conds",
                                                 freqs.ending="freqs" ){

    # read in or set default background frequencies
    if( useFreqs ){
        file.freqs <- paste( filename, freqs.ending, sep="." )
        freqs <- scan( file.freqs, quiet=TRUE )
    } else{
        freqs <- rep( .25, 4 )
    }

    # read in conditional probabilities
    file.conds <- paste( filename, conds.ending, sep="." )
    conds <- readBaMM( file.conds )

    # read in probabilities
    file.probs <- paste( filename, probs.ending, sep="." )
    probs <- readBaMM( file.probs )

    MAXORDER <- 5
    
    # model formats
    L <- unique( c( sapply( probs, nrow ), sapply( conds, nrow ) ) )
    maxorder <- unique( length( probs ), length( conds ) ) - 1
    if( length( L ) > 1 || length( maxorder ) > 1 ){
        stop( paste( file.probs, " and ", file.conds, " formats differ.", sep=""
              ) )
    }
    
    # determine maximum order
    maxorder <- min( maxorder, MAXORDER )

    # calculate the position-specific contributions of k-mers to the information
    # content
    kmerContributions <- kmerContributionsToInformationContent( probs, conds, 0:maxorder, freqs=freqs, base=base, rna=rna )

    # calculate the position-specific information content from k-mer
    # contributions
    informationContent <- informationContent.kmerContributions( kmerContributions, unsigned=TRUE )

    # create color palette
    diagram.color<- rgb( matrix( c( c( 0, 69, 134 ), c( 255, 66, 14 ),
                         c( 255, 211, 32 ), c( 87, 157, 28 ), c( 126, 0, 33 ),
                         c( 131, 202, 255 ), c( 49, 64, 4 ), c( 174, 207, 0 ),
                         c( 75, 31, 111 ), c( 255, 149, 14 ), c( 197, 0, 11 ),
                         c( 0, 132, 209 ) ), ncol=3, byrow=TRUE ), max=255 )
    # sample colors
    color <- diagram.color[ sample( 1:length( diagram.color ), maxorder+1 ) ]

    # plot higher-order contributions
    barplot( informationContent, names.arg=1:L, legend.text=0:maxorder,
             beside=TRUE, col=color, border=color, ylim=c( 0, 2 ),
             xlab="Model position [nt]", ylab="Information content",
             args.legend=list( x="topright", border=color, box.col=par( "bg" ),
             inset=.02, title="Order" ) )
}

totalContributionsToInformationContent <- function( filename, order,
                                                    useFreqs=FALSE, base=2,
                                                    rna=FALSE,
                                                    probs.ending="probs",
                                                    conds.ending="conds",
                                                    freqs.ending="freqs" ){

    # read in or set default background frequencies
    if( useFreqs ){
        file.freqs <- paste( filename, freqs.ending, sep="." )
        freqs <- scan( file.freqs, quiet=TRUE )
    } else{
        freqs <- rep( .25, 4 )
    }

    # read in conditional probabilities
    file.conds <- paste( filename, conds.ending, sep="." )
    conds <- readBaMM( file.conds )

    # read in probabilities
    file.probs <- paste( filename, probs.ending, sep="." )
    probs <- readBaMM( file.probs )

    MAXORDER <- 5

    # model formats
    L <- unique( c( sapply( probs, nrow ), sapply( conds, nrow ) ) )
    maxorder <- unique( length( probs ), length( conds ) ) - 1
    if( length( L ) > 1 || length( maxorder ) > 1 ){
        stop( paste( file.probs, " and ", file.conds, " formats differ", sep=""
              ) )
    }
    
    # maximum model order
    if( any( order > maxorder ) ){
        missing.order <- order[order > maxorder]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " missing", sep="" ) )
    }

    # maximum code order
    if( any( order > MAXORDER ) ){
        missing.order <- order[order > MAXORDER]
        p <- length( missing.order ) > 1
        warning( paste( "Order", ifelse( p, paste( "s ", paste( rev( rev(
                 missing.order )[-1] ), collapse=", " ), " and ", rev(
                 missing.order )[1], sep="" ), paste( " ", missing.order, sep=""
                 ) ), " not implemented yet", sep="" ) )
    }

    # determine maximum order
    maxorder <- min( maxorder, MAXORDER, max( order ) )

    # k-mer contributions
    kmerContr <- kmerContributionsToInformationContent( probs, conds, order,
                                                        freqs, base, rna )
    # positional contributions
    posContr <- informationContent.kmerContributions( kmerContr, unsigned=TRUE )
    
    # order contributions
    orderContr <- apply( posContr, 1, sum, na.rm=TRUE )

    return( orderContr )
}
