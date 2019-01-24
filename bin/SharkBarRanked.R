library("reshape2")
library("ggplot2")

dd <- read.table('SharkBayVHighbourneForBar.txt', header = 1, sep="\t")
mm <- melt(dd)

# refactor "variable" aka names by the order I want:
order_i_want <- mm[with(mm, order(-value)),]$variable
mm$variable = factor(mm$variable, levels=rev( unique(order_i_want)) ) 

#Plot ggplot2 with coord_flip
p2 <- ggplot(mm,aes(x=factor(variable), y=value, fill=factor(Taxonomy))) + geom_bar(position=position_dodge()) + coord_flip() + xlab("Taxa") + ylab("Relative Abundance (%)") + guides(fill=guide_legend(title="")) + theme_bw()
