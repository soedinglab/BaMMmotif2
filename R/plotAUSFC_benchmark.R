#!/usr/bin/Rscript
#-----------------------------------------------------------------
#
# Documentation
#
#-----------------------------------------------------------------

# read in FDR and Sensitivity values from .zoops.stats file 
# plot FDR vs. Sensitivity curve
# caluclate area under the sensitivity-FDR curve(AUSFC)
# AUSFC can be used for scoring large datasets for benchmarking
# AUSFC is saved in a .ausfc file

#-----------------------------
#
# R sources
#
#-----------------------------

# load "zoo" package for calculating AUSFC
if (!require("zoo")){
  install.packages("zoo", repos="http://cran.rstudio.com/")
}

#-----------------------------
#
# Parameter setting
#
#-----------------------------

args <- commandArgs(trailingOnly=TRUE)

# get the directory to p-value files
dir <- args[1]
dataname <- args[2]

# local test:
#dir = c("/home/wanwan/benchmark/JunD/EM_FDR_benchmark/")
#dataname = c("JunD_motif_1")

# read in the .stats files
stat <- read.table(paste(dir, dataname, ".zoops.stats", sep = "" ), skip=1)
fdr <- stat$V3 
recall <- stat$V4

len = length(fdr)
list <- seq(1, len)

# avoid FDR from being zero
for(i in list){
  if( fdr[i] ==0 ){
    fdr[i] = 0.000001
  }
}

# modify the SF curve:
# reset recall to 1 when it is larger than 1
for(i in list){
  if( recall[i] > 1 ){
    for(rest in i:len){
      recall[rest] = 1
    }
    break
  }
}

# extend the SF curve till FDR reaches 0.5
recall <- append(recall, 1)
fdr <- append(fdr, 0.5)
list <- seq(1, len+1)

# limit the frame of the curve
left = 1
right = len+1

for(i in list){
  if( fdr[i] >= 0.00009 ){
    left = i
    break
  }
}

for(i in list){
  if(fdr[i] >= 0.5 ){
    right = i
    break
  }
}

list <- seq(left, right)

# plot fdr vs. recall curve
sum_area = log10(0.5) + 4
png( file = paste( dir, dataname, '_SFCurve_cpp.png', sep = "" ) )
plot(fdr[list], recall[list], log="x", 
     main=paste(dataname, "Sensitivity vs. FDR", sep="\n\n"),
     xlab="FDR", ylab="Sensitivity", xlim=c(0.0001,0.5), ylim=c(0,1),
     type='l', lwd=2.5,
     col="deepskyblue")
# compute  the area under the sensitivity-FDR curve(AUSFC):
ausfc = sum(diff(log10(fdr[list]))*rollmean(recall[list],2)) / sum_area
# write the AUSFC on the plot with the precision of 3 digits:
text(0.0005, 0.95, paste( "AUSFC: ", round(ausfc, digits=4) ))
dev.off()

# write paramaters to a file
ofile = paste( dir, dataname, '.ausfc', sep = "" )
writeLines( c( "AUSFC: ", ausfc ), ofile )
