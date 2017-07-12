#!/usr/bin/env Rscript
#####################################################
# Plot three curves for false discovery rate estimation
# 1. The (TP+FP)+ versus FP- curve for the MOPS
# 2. The precision versus recall curve for the MOPS model
# 3. The precision versus recall curve for the ZOOPS model
# 4. Comparision of precision-recall curves for bothe MOPS and ZOOPS models
#####################################################

args <- commandArgs( trailingOnly = TRUE )
print( args )

dir <- args[1]		# Directory for input and output
name <- args[2]		# name of the dataset 

# plot FP vs. TP+FP curve for MOPS model
FPmops <- read.table( paste( dir, name, ".mops.fp", sep = "" ), fileEncoding="latin1" )
TFPmops <- read.table( paste( dir, name, ".mops.tfp", sep = "" ), fileEncoding="latin1" )
png( file = paste( dir, name, '_TFP_vs_FP_MOPS.png', sep = "" ) )
par( oma = c( 1,2,1,1 ) )
plot(
	FPmops$V1, TFPmops$V1, 
	main = paste( name, "(MOPS model)" ), 
	xlab = "FP-", 
	ylab = "(TP+FP)+", 
	type = "l", lwd = 2.5,
	col = "chocolate1",
	cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5
)
abline( 0, 1, lwd = 2.5, col = "grey66" )
dev.off()

# plot precision vs. recall curve for MOPS model
recallMOPS <- read.table( paste( dir, name, ".mops.recall", sep = "" ), fileEncoding="latin1" )
precisionMOPS <- read.table( paste( dir, name, ".mops.precision", sep = "" ), fileEncoding="latin1" )
png( file = paste( dir, name, '_precision_recall_MOPS.png', sep = "" ) )
par( oma = c(1,1,1,1) )
plot(
	recallMOPS$V1, precisionMOPS$V1, 
	main = paste( name, "(MOPS model)" ), 
	#xlim = range( 0,1 ), ylim = range( 0,1 ), 
	xlab = "recall (sensitivity)", 
	ylab = "precision (positive predictive value)", 
	type = "l", lwd = 2.5,
	col = "deepskyblue",
	cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5
)
dev.off()

# plot precision vs. recall curve for ZOOPS model
recallZOOPS <- read.table( paste( dir, name, ".zoops.recall", sep = "" ), fileEncoding="latin1" )
precisionZOOPS <- read.table( paste( dir, name, ".zoops.precision", sep = "" ), fileEncoding="latin1" )
png( file = paste( dir, name, '_precision_recall_ZOOPS.png', sep = "" ) )
par( oma=c(1,1,1,1) )
plot(
	recallZOOPS$V1, precisionZOOPS$V1, 
	main = paste( name, "(ZOOPS model)" ), 
	xlim = range(0,1), ylim=range(0,1), 
	xlab = "recall (sensitivity)", 
	ylab = "precision (positive predictive value)", 
	type = "l", lwd = 2.5,
	col = "chartreuse4",
	cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5
)
dev.off()

# Merge precision-recall curves for both MOPS and ZOOPS models on one plot
png( file = paste( dir, name, '_precision_recall_comparison.png', sep = "" ) )
par( oma = c(1,1,1,1), mar = c(5,5,7,5) )
### Plot for MOPS:
plot(
	recallMOPS$V1, precisionMOPS$V1, 
	main = paste( "precision-recall curves\n for MOPS and ZOOPS models \n(", name, ")" ), 
	xlim = range(0,1), ylim = range(0,1), 
	xlab = "recall (sensitivity)", 
	ylab = "precision (positive predictive value)", 
	type = "l", lwd = 2.5,
	col = "deepskyblue",
	cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, cex.sub = 1.5
)
par( new = T )
#### Plot for ZOOPS: 
plot(
	recallZOOPS$V1, precisionZOOPS$V1, 
	xlim = range(0,1), ylim = range(0,1), 
	xlab = "", 
	ylab = "", 
	type = "l", lwd = 2.5,
	col = "chartreuse4",
	axes = F
)
par( new = F )
legend(
	'bottomleft', 
	c( "MOPS model", "ZOOPS model" ), 
	lty = 1, lwd = 2.5, 
	col = c( "deepskyblue", "chartreuse4" ), 
	bty = 'n', cex = 1
)
dev.off()

