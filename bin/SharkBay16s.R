#Dot Plot - SharkBay 16S
# Written RAWIII April 12, 2014

#Set Working Directory
setwd("/home/rawiii/Desktop/SharkBayR/")

#Load Libraries
library("ggplot2")
library("RColorBrewer")
library("reshape")
library("reshape2")
library("lattice")

#Load Data and Rename Columns
data.16s <- read.table("SB16SBacformat.txt", sep="\t", header=T)
data.16s.m <- melt(data.16s)
data.18s <- read.table("SB18SEukformat.txt", sep="\t", header=T)
data.18s.m <- melt(data.18s)
data.refseq.bac <- read.table("SBrefseqBacformat.txt", sep="\t", header=T)
data.refseq.bac.m <- melt(data.refseq.bac)
data.refseq.Euk <- read.table("SBrefseqEukformat.txt", sep="\t", header=T)
data.refseq.Euk.m <- melt(data.refseq.Euk)

#Plot Dot Plot
svg('/home/rawiii/Desktop/SharkBayR/Test16s.svg', width = 18, height = 12, pointsize = 16)
p1 <- ggplot(data.16s.m, aes(x=variable, y=phylum, color=variable)) + geom_point(aes(size=log(value))) +  theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p2 <- ggplot(data.refseq.bac.m, aes(x=variable, y=phylum, color=variable)) + geom_point(aes(size=log(value))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
grid.arrange(arrangeGrob(p1,p2, ncol=3, widths=c(5/7, 5/7, 1/16)))
dev.off()


svg('/home/rawiii/Desktop/SharkBayR/Test18s.svg', width = 18, height = 18, pointsize = 16)
p3 <- ggplot(data.18s.m, aes(x=variable, y=phylum, color=variable)) + geom_point(aes(size=log(value))) +  theme_bw() + theme(axis.text.x = element_text())
p4 <- ggplot(data.refseq.Euk.m, aes(x=variable, y=phylum, color=variable)) + geom_point(aes(size=log(value))) + theme_bw() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
grid.arrange(arrangeGrob(p3,p4, ncol=1, widths=c(5/7, 5/7, 1/16)))
dev.off()

#Cluster Euclidian
data.t  <- read.table("SB16SformatEdit.txt", sep="\t", header=T, row.names=1)
data.tt  <- t(data.t)
data.tt.dist <- dist(data.tt)
data.tt.euclid.fit  <- hclust(data.tt.dist)
plot(data.tt.euclid.fit, main="Distance: Euclidian, Clustering: Complete")

#pvclust with monte carlo euclidian distance
library(pvclust)
data.pv_fit <- pvclust(data.t, method.hclust="complete", method.dist="euclidian", n=1000)
plot(data.pv_fit)

#MDS plot
library(vegan)
data.t  <- read.table("SB16SformatEdit.txt", sep="\t", header=T, row.names=1)
data.hel  <- sqrt(data.t)
data.hel.pca  <- rda(t(data.hel))
p <- length(data.hel.pca$CA$eig)
data.hel.pca.sc1 <- scores(data.hel.pca, display="wa", scaling=1, choices=c(1:p))
variance = (data.hel.pca$CA$eig / sum(data.hel.pca$CA$eig))*100
plot(data.hel.pca.sc1, type="p", xlab= paste("PC1 (", round(variance[1],2) ,"% Variance)"), 
     ylab= paste("PC2 (", round(variance[2],2) ,"% Variance)"),)

#Plot MDS ggplot2
quartz("Pathways Scaling 1: PCA")
qplot(data.hel.pca.sc1[,1], 
      data.hel.pca.sc1[,2], 
      label=rownames(data.hel.pca.sc1), 
      size=2, geom=c("point"), 
      xlab= paste("PC1 (", round(variance[1],2) ," % Variance)"), 
      ylab= paste("PC2 (", round(variance[2],2) ," % Variance)"), 
      color=factor(data.hel.pca)) + 
  geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() + theme(legend.position="none") + xlim(-0.6, 0.6)

#MDS plot
library(vegan)
data.hel.nmds  <- metaMDS(t(data.hel), distance = "bray")
qplot(data.hel.nmds$points[,1], data.hel.nmds$points[,2], label=rownames(data.hel.nmds$points), size=2, geom=c("point"), 
      xlab="MDS1", ylab="MDS2", main=paste("NMDS/Bray - Stress =", round(pdata.hel.nmds$stress,3)), color=factor(data.hel.nmds)) + 
  geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() +theme(legend.position="none") + xlim(-0.5,1.0)