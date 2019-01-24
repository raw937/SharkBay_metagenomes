#PCA refseq SharkBay May 3, 2013

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
data.t  <- read.table("SBrefSeqformat.txt", sep="\t", header=T, row.names=1)
data.hel  <- sqrt(data.t)
data.hel.pca  <- rda(t(data.hel))
p <- length(data.hel.pca$CA$eig)
data.hel.pca.sc1 <- scores(data.hel.pca, display="wa", scaling=1, choices=c(1:p))
variance = (data.hel.pca$CA$eig / sum(data.hel.pca$CA$eig))*100

#cluster groups
source("pvclust_bcdist.R")
pathways_wide.hel.pv_fit <- pvclust(as.matrix(data.t), method.hclust="ward", method.dist="bray–curtis", n=1000)

# look at fit and decide cut height
plot(pathways_wide.hel.pv_fit) 
pathways_wide.hel.groups <- cutree(pathways_wide.hel.pv_fit$hclust, h=1.25) # slice dendrogram for groups

#Plot PCA
library(grid)
library(gridExtra)
p1 <- qplot(data.hel.pca.sc1[,1], 
            data.hel.pca.sc1[,2], 
            label=rownames(data.hel.pca.sc1), 
            size=2, geom=c("point"), 
            xlab= paste("PC1 (", round(variance[1],2) ," % Variance)"), 
            ylab= paste("PC2 (", round(variance[2],2) ," % Variance)"), 
            color= factor(pathways_wide.hel.groups)) + 
  geom_text(hjust=-0.1, vjust=0, colour="black", size=3) + theme_bw() + xlim(-11,7) + theme(legend.position="none") 

p2  <- qplot(data.hel.pca.sc1[,1], 
             data.hel.pca.sc1[,3], 
             label=rownames(data.hel.pca.sc1), 
             size=2, geom=c("point"), 
             xlab= paste("PC1 (", round(variance[1],2) ," % Variance)"), 
             ylab= paste("PC3 (", round(variance[3],3) ," % Variance)"), 
             color= factor(pathways_wide.hel.groups)) + 
  geom_text(hjust=-0.02, vjust=0, colour="black", size=3) +  theme_bw() + xlim(-11,7) + theme(legend.position="Test")

p3  <- qplot(data.hel.pca.sc1[,2], 
             data.hel.pca.sc1[,3], 
             label=rownames(data.hel.pca.sc1), 
             size=2, geom=c("point"), 
             xlab= paste("PC2 (", round(variance[2],2) ," % Variance)"), 
             ylab= paste("PC3 (", round(variance[3],3) ," % Variance)"), 
             color= factor(pathways_wide.hel.groups)) + 
  geom_text(hjust=-0.02, vjust=0, colour="black", size=3) + theme_bw() + xlim(-11,7) + theme(legend.position="Test")

grid.arrange(p1, p2, p3, ncol=1, main = "")

