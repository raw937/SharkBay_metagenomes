# indicator species analysis
library('indicspecies')


data.t  <- read.table("SB16SformatEdit.txt", sep="\t", header=T, row.names=1)
ind_val_data <- t(data.t) # samples as rows
groups = c(rep(1, 3),2 , rep(3,3), 2,2) # rows (samples) that are related
indval <- multipatt(as.data.frame(ind_val_data), groups, control = how(nperm=999))
summary(indval) # basic overview

# sign gives full table of indicator species
indval$sign

# plot p-value and decide on cutoff
plot(sort(indval$sign$p.value))
# i.e. p<0.1
indval$sign[indval$sign$p.value < 0.1,]