#Load Libraries
library("vegan")
library("ggplot2")
library("RColorBrewer")
library("reshape")
library("reshape2")
library("lattice")
library("ecodist")
library("rgl")


#Load data
data.t  <- read.table("SB16SformatEdit.txt", sep="\t", header=T, row.names=1)
data.hel  <- sqrt(data.t)
data.hel.pca  <- rda(t(data.hel))
p <- length(data.hel.pca$CA$eig)
data.hel.pca.sc1 <- scores(data.hel.pca, display="wa", scaling=1, choices=c(1:p))
variance = (data.hel.pca$CA$eig / sum(data.hel.pca$CA$eig))*100

# cluster groups
try(library("devtools"), install.packages("devtools")) # used to source functions from the internet
library("devtools")
source_url('http://raw.github.com/nielshanson/mp_tutorial/master/taxonomic_analysis/code/pvclust_bcdist.R')

source("pvclust_bcdist.R")
pathways_wide.hel.pv_fit <- pvclust(as.matrix(data.t), method.hclust="ward", method.dist="bray–curtis", n=1)
plot(pathways_wide.hel.pv_fit) # look at fit and decide cut height
pathways_wide.hel.groups <- cutree(pathways_wide.hel.pv_fit$hclust, h=1.0) # slice dendrogram for groups

#Plot MDS
qplot(data.hel.pca.sc1[,1], 
      data.hel.pca.sc1[,2], 
      label=rownames(data.hel.pca.sc1), 
      size=2, geom=c("point"), 
      xlab= paste("PC1 (", round(variance[1],2) ," % Variance)"), 
      ylab= paste("PC2 (", round(variance[2],2) ," % Variance)"), 
      color= factor(pathways_wide.hel.groups)) + 
  geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() + xlim(-2,2) + theme(legend.position="none") 

# what are the primary loadings of each PC?
sort(abs(data.hel.pca$CA$v.eig[,1]), decreasing=TRUE)
sort(abs(data.hel.pca$CA$v.eig[,2]), decreasing=TRUE)
sort(abs(data.hel.pca$CA$v.eig[,3]), decreasing=TRUE)

# plot the e.vectors plot loadings for the first three PCs
plot(sort(abs(data.hel.pca$CA$v.eig[,1]), decreasing=TRUE))
points(sort(abs(data.hel.pca$CA$v.eig[,2]), decreasing=TRUE), col="red")
points(sort(abs(data.hel.pca$CA$v.eig[,3]), decreasing=TRUE), col="blue")

# being informed from the loadings, one might look at proteos and cyanos directly
plot(t(data.t[c("Proteobacteria", "Cyanobacteria"),]), col=factor(pathways_wide.hel.groups))

# use the pairs to plot the guys different
pairs(t(data.t[c("Proteobacteria", "Unclassified Eukaryota", "Arthropoda", "Cyanobacteria"),]), col=factor(pathways_wide.hel.groups))


plot(variance) # look for the elbow in scree plot

#plot 3D PCA
library(rgl)
plot3d(data.hel.pca.sc1[,1], 
       data.hel.pca.sc1[,2], 
      data.hel.pca.sc1[,3], col=factor(pathways_wide.hel.groups), size=10)


#MDS plot
library(vegan)
data.hel.nmds  <- metaMDS(t(data.hel), distance = "bray")
qplot(data.hel.nmds$points[,1], data.hel.nmds$points[,2], label=rownames(data.hel.nmds$points), size=2, geom=c("point"), 
      xlab="MDS1", ylab="MDS2", main=paste("NMDS/Bray - Stress =", round(data.hel.nmds$stress,3)), color=factor(pathways_wide.hel.groups)) + 
  geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() +theme(legend.position="none") 