#GC plot from Perl Script Output - SharkBay Contigs
# Written RAWIII Oct 2, 2013

#Set Working Directory
setwd("/home/rawiii/Desktop/SharkBayR/")

#Load Libraries
library("ggplot2")
library("RColorBrewer")
library("reshape")
library("reshape2")
library("lattice")

#Load Data and Rename Columns
data <- read.table("SharkBay_GCcount.txt")
names(data) = c("Sample", "Contig_ID", "GC")

#Plot with Density ggPlot2
m <- ggplot(data, aes(x=GC, colour=Sample, group=Sample))
m = m + geom_density(fill=NA)
m
