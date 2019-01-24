# Shark Bay MetaPathways MetaCyc Venn Diagram and ggplot2 Pathway Table
# Written by RAWIII Sep 12 2013

# Load Libraries
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("reshape")
library("reshape2")

# Set working directory
setwd('/home/rawiii/Desktop/SharkBayR/')

# Source for R program operation
source("venn_diagram3_wSvg.r")

# Load MetaCyc v17 database
meta17 <- read.table("meta_17.txt", sep="\t", header=T, row.names=1)

# Load MetaCyc v17 Hierarchy and set levels
meta_17_hier <- read.table("meta_17_hierarchy.txt", sep="\t", row.names=1)
meta_17_hier$V3 <- factor(meta_17_hier$V3, levels=union(meta_17_hier$V3,meta_17_hier$V3))
meta_17_hier$V4 <- factor(meta_17_hier$V4, levels=union(meta_17_hier$V4,meta_17_hier$V4))

# Load Datasets (Delim)
Columnar  <- read.delim("Columnar.txt", row.names=1)
Pustular  <- read.delim("Pustular.txt", row.names=1)
Smooth  <- read.delim('Smooth.txt', row.names=1)

# Replace with Simpler header
header <- c("long_name", "num_rxn", "num_covered", "num_orfs")
colnames(Smooth) <- header
colnames(Pustular) <- header
colnames(Columnar) <- header

# Trim pathways that do not have 5 Orfs
Columnar <- droplevels(Columnar[Columnar["num_orfs"] >5,])
Pustular <- droplevels(Pustular[Pustular["num_orfs"] >5,])
Smooth <- droplevels(Smooth[Smooth["num_orfs"] >5,])

normalize_length = TRUE # normalize by pathway length
normalize_orfs = TRUE # normalize to number of ORFs in sample

# total Orfs predicted for each sample 
num_orfs_Col = 54839
num_orfs_Pus = 248416
num_orfs_Smo = 111655

# normalize all pathway counts to length
if (normalize_length == TRUE) {
Columnar[,"num_orfs"] <- Columnar[,"num_orfs"] / Columnar[,"num_rxn"]
Pustular[,"num_orfs"] <- Pustular[,"num_orfs"] / Pustular[,"num_rxn"]
Smooth[,"num_orfs"] <- Smooth[,"num_orfs"] / Smooth[,"num_rxn"]
}

if (normalize_orfs == TRUE) {
Columnar[,"num_orfs"] <- (Columnar[,"num_orfs"] / num_orfs_Col) * 100
Pustular[,"num_orfs"] <- (Pustular[,"num_orfs"] / num_orfs_Pus) * 100
Smooth[,"num_orfs"] <- (Smooth[,"num_orfs"] / num_orfs_Smo) * 100
}

# Plot Venn Diagram 
SharkBay_venn3  <- venn_diagram3(rownames(Columnar), rownames(Pustular), rownames (Smooth), "Columnar", "Pustular", "Smooth", name_output = "SharkBayVenn")

# Construct a dataframe for Top_50
Columnar_50 <- Columnar[ order(Columnar$num_orfs, decreasing=TRUE), ][1:50,]
Pustular_50 <- Pustular[ order(Pustular$num_orfs, decreasing=TRUE), ][1:50,]
Smooth_50 <- Smooth[ order(Smooth$num_orfs, decreasing=TRUE), ][1:50,]

# Create union of all pathways
union_top50 <- union( union( rownames(Columnar_50), rownames(Pustular_50) ), rownames(Smooth_50))

#Isolate pathways in data.frame object
top_50 <- data.frame( pathway=union_top50, long_name=meta17[union_top50,], hier1 = meta_17_hier[union_top50, 2], hier2 = meta_17_hier[union_top50, 3], Columnar=rep(0,length(union_top50)), Pustular=rep(0,length(union_top50)), Smooth=rep(0,length(union_top50)) )

#Fill in data
rownames(top_50) = union_top50
top_50[union_top50,"Columnar"] = Columnar[union_top50,"num_orfs"]
top_50[union_top50,"Pustular"] = Pustular[union_top50,"num_orfs"]
top_50[union_top50,"Smooth"] = Smooth[union_top50,"num_orfs"]

#Fill in missing data with zeros
top_50[is.na(top_50)] <- 0

#Pathway Order based on Hier1 and Hier2 Hierchy Levels
pathway_ordering = rev(top_50[ with(top_50,order(hier1,hier2)), "long_name"])
top_50$long_name <- factor(top_50$long_name, levels=pathway_ordering)

#Use melt function to prepare for ggplot in long data form
top_50_m <- melt(top_50, id.vars = c("pathway","long_name"), measure.vars = c("Columnar", "Pustular", "Smooth"))
top_50_m[is.na(top_50_m)] <- 0

#Drop extra levels if any
top_50_m <- droplevels(top_50_m)

#Perform Plot
ggplot(top_50_m, aes(long_name, value, fill=variable)) + geom_bar() + facet_wrap(~ variable, ncol=4) + opts(axis.text.x=theme_text(angle=90)) +
  scale_fill_manual(values=c(rgb(90,184,75, max = 255),rgb(113,200,199, max = 255),rgb(220,86,232, max = 255),rgb(210,7,62, max = 255))) + coord_flip()
