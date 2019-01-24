#Dot Plot - RefSeq PL surr
# Written RAWIII May 5, 2014

#Set Working Directory
setwd("/home/rawiii/Desktop/SharkBayR/")

#Load Libraries
library("ggplot2")
library("RColorBrewer")
library("reshape")
library("reshape2")
library("lattice")

#Load Data and Rename Columns
data <- read.table("PL_surr_refseq_top25norm.txt", sep="\t", header=T)
data.m <- melt(data)

ggplot(data.m, aes(x=variable, y=Phyla, color=variable)) + geom_point(aes(size=log(value))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
