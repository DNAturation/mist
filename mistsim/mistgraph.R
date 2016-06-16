#draws graph for mistsim csv file

library(optparse)
library(ggplot2)
library(gridExtra)

args<-commandArgs(trailingOnly=T)
#make_option('-d', '--directory', default='./sim/
#make_option('-i', '--simfile', default='simulator.csv', help='.csv file made by mistsim')
#make_option('-o', '--outfile', default='Simgraph.png', help='name of output file')

setwd(args[1])

graphdata <- read.csv(args[2], header=TRUE, sep=" ")


plot1<-ggplot(graphdata, aes(x = GenomeThreshold, y=GenomesLeft)) + geom_line() +xlab('Less than X number of
genes missing from genomes') +ggtitle('plot of genomes remaining with threshold cutoff')
plot2<-ggplot(graphdata, aes(x = GeneThreshold, y=GenesLeft)) + geom_line() +xlab('Less than X number of
genomes missing from genes') + ggtitle('plot of genes remaining with threshold cutoff')
png(args[3])
grid.arrange(plot1, plot2)

garbage<-dev.off()

