################################################################################
#
# Documentation
# 
################################################################################

# Plot the (positional) k-mer(s) that contribute(s) most to the information
# content
# using files with (homogeneous) BaMM probabilities as obtained from BaMM!motif

################################################################################
#
# R Sources
# 
################################################################################

source( "plotBaMM.R" )
source( "readBaMM.R" )

################################################################################
#
# Code to plot max. contributing k-mer(s)
# 
################################################################################

plotMaxContributingKmer <- function( filename, order, base=2, rna=FALSE,
                                     useFreqs=FALSE, probs.ending="probs",
                                     conds.ending="conds", freqs.ending="freqs"
                                     ){

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
    
    # determine order
    order <- min( order, maxorder, MAXORDER )

    # calculate the position-specific contributions of k-mers to the information
    # content
    kmerContributions <- kmerContributionsToInformationContent( probs, conds,
                             order, freqs=freqs, base=base, rna=rna )[[order+1]]

    # calculate the max. position-specific k-mer contribution to the information
    # content
    maxKmerContributions <- apply( kmerContributions[,-( 1:( order+1 ) ),
                                   drop=FALSE], 2, max )

    ylim <- c( min( kmerContributions[,-( 1:( order+1 ) ), drop=FALSE] ), max(
               maxKmerContributions ) + .3 * max( maxKmerContributions ) )

    linesplot( kmerContributions[,-( 1:( order+1 ) ), drop=FALSE], ylim=ylim,
               main="", xlab="Position", ylab="Information Content" )

    # extract the max.-contributing k-mers
    kmer <- as.character( apply( kmerContributions[apply( kmerContributions[,-(
                          1:( order+1 ) ), drop=FALSE], 2, which.max ), 1:(
                          order+1 ), drop=FALSE], 1, paste, collapse="" ) )

    bases <- c( "A", "C", "G", ifelse( rna, "U", "T" ) )
    indices <- 1:4
    names( indices ) <- bases

    color.bases <- c( "#008000", "#0000FF", "#FFA500", "#FF0000" )

    for( i in 1:( order+1 ) ){
        a <- substr( kmer, 1, i-1 )
        b <- substr( kmer, i, i )
        c <- substr( kmer, i+1, order+1 )
        labels <- apply( cbind( a, b, c ), 1, function( x ){ bquote( phantom( .(
                         x[1] ) ) * .( x[2] ) * phantom( .( x[3] ) ) ) } )
        colors <- color.bases[ as.integer( indices[ b ] ) ]
        text( ( 1:( L-order ) )+.1, pos=3, maxKmerContributions, labels=sapply(
              labels, as.expression ), offset=2, col=colors, srt=90 )
    }
}

plotMaxContributingBgKmer <- function( filename, order, base=2, rna=FALSE,
                                       probs.ending="probs",
                                       conds.ending="conds" ){

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
        stop( paste( file.probs, " and ", file.conds, " formats differ.", sep=""
              ) )
    }

    # determine order
    order <- min( order, maxorder, MAXORDER )

    # calculate the contributions of k-mers to the information content
    kmerContributions <- bgKmerContributionsToInformationContent( probs, conds,
                             order, freqs=freqs, base=base, rna=rna )[[order+1]]

    # calculate the max. k-mer contribution to the information content
    maxKmerContributions <- max( kmerContributions[,-( 1:( order+1 ) )] )

    ylim <- c( min( kmerContributions[,-( 1:( order+1 ) )] ), max(
               maxKmerContributions ) + .3 * max( maxKmerContributions ) )

    linesplot( kmerContributions[,-( 1:( order+1 ) )], labels=order, ylim=ylim,
               main="", xlab="Order", ylab="Information Content" )

    # extract the max. contributing k-mer
    kmer <- as.character( paste( kmerContributions[which.max( kmerContributions[,-(
                          1:( order+1 ) )] ), 1:( order+1 ) ], collapse="" ) )

    bases <- c( "A", "C", "G", ifelse( rna, "U", "T" ) )
    indices <- 1:4
    names( indices ) <- bases

    color.bases <- c( "#008000", "#0000FF", "#FFA500", "#FF0000" )

    for( i in 1:( order+1 ) ){
        a <- substr( kmer, 1, i-1 )
        b <- substr( kmer, i, i )
        c <- substr( kmer, i+1, order+1 )
        labels <- apply( cbind( a, b, c ), 1, function( x ){ bquote( phantom( .(
                         x[1] ) ) * .( x[2] ) * phantom( .( x[3] ) ) ) } )
        colors <- color.bases[ as.integer( indices[ b ] ) ]
        text( ( 1:( L-order ) ), pos=3, maxKmerContributions, labels=sapply(
              labels, as.expression ), offset=2, col=colors, srt=90 )
    }
}

################################################################################
#
# This code is built on code from the CRAN package LSD
# https://cran.r-project.org/web/packages/LSD/index.html
# 
################################################################################

linesplot <- function( x, labels=NULL, col="black", cols=NULL, alpha=25,
                       ylim=NULL, xlab=NULL, ylab="data values", las=1,
                       outline=TRUE, cexbox=0.6, addboxes=FALSE, border="black",
                       range=1.5, lwd=1.5, main="Linesplot", ... ){

    if( !( is.vector( x ) ) & !( is.matrix( x ) ) & !( is.list( x ) ) & !(
        is.data.frame( x ) ) ){
        stop( "x must be a vector, matrix, list or a data.frame" )
    }
    if( !( is.list( x ) ) & !( is.matrix( x ) ) & !( is.data.frame( x ) ) ){
        x <- cbind( x )
    }
    if( is.data.frame( x ) ){
        x <- as.list( x )
    }
    if( is.matrix( x ) ){
        x <- as.list( as.data.frame( x ) )
    }
    if( is.null( labels ) ){
        labels <- labels( x )
    }
    if( is.null( cols ) ){
        cols <- rep( convertcolor( col, alpha ), length( x ) )
    } else{
        cols <- convertcolor( cols, alpha )
    }
    if( is.null( ylim ) ){
        ylim <- range( unlist( x ), na.rm=TRUE )
    }
    if( length( x ) == 1 ){
        cexbox <- 0.4
    }
    if( length( x ) > 1 ){
        boxwex <- cexbox
    } else{
        boxwex <- NULL
    }
    plot( 1, col="white", xlim=c( 0.5, length( x ) + 0.5 ), ylim=ylim,  main="",
          xlab=xlab, ylab=ylab, xaxt="n", ... )
    par.axis.default <- par( "cex.axis" )
    par( cex.axis=par( "cex.lab" ) )
    axis( 1, at=seq( 1, length( x ), 1 ), labels=labels, las=las )
    for( j in 1:length( x ) ){
        for( i in x[[j]] ){
            lines( c( j-cexbox/2, j+cexbox/2 ), c( i, i ), lwd=2, col=cols[j] )
        }
    }
    if( addboxes ){
        boxplot( x, add=TRUE, col="transparent", xlim=c( 0.5, length( x )+0.5 ),
                 ylim=ylim, width=NULL, boxwex=boxwex, outline=outline,
                 axes=FALSE, border=border, lwd=lwd, range=range )
    }
    par( cex.axis = par.axis.default )
}
