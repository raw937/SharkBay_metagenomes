library("ggplot2")
library("reshape2")

d <- read.delim("SharkvHighbourneViromes.txt")

colnames(d) <- c("Taxa", "Shark Bay","Highbourne Cay")
d <- melt(d)

ggplot(d, aes(Taxa)) + geom_bar() + facet_wrap(~ value)

ggplot(d, aes(Taxa, fill=variable)) + geom_bar(position="dodge")

p1<-ggplot(mm,aes(x=factor(Taxa),y=value,fill=factor(variable)), color=factor(variable)) + geom_bar(position=position_dodge())



mm <- d

# refactor "variable" aka names by the order I want:
order_i_want <- mm[with(mm, order(-value)),]$variable
mm$variable = factor(mm$variable, levels=rev( unique(order_i_want)) )

#Plot ggplot2 with coord_flip
ggplot(mm,aes(x=factor(variable), y=value, fill=factor(Taxa))) + geom_bar(position=position_dodge()) + coord_flip() + xlab("SEED Subsystem") + ylab("Number of Hits") + guides(fill=guide_legend(title=""))
