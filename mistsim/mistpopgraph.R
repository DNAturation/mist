#draws graph for mistpopper

library(optparse)
library(ggplot2)
library(gridExtra)
args<-commandArgs(trailingOnly=T)
setwd(args[1])
graphdata <- read.csv(args[2], header=TRUE, sep=" ")

png(args[3])

ggplot(graphdata, aes(x = NumberOfGenesRemoved, y=NumberOfPerfectGenomes)) + geom_line(
) +xlab('number of genes popped') +ylab('number of perfect genomes')+ggtitle('plot of perfect genomes after removing worst X number of genes')

garbage<-dev.off()