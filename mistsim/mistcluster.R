library(ape)
library(data.table)
library(tools)
library(gtools)
library(gridExtra)
source('wallace-helper.R')
args<-commandArgs(trailingOnly=T)

#1 = tables? .matrix extension
#2 = outfile?
cutnames<-file_path_sans_ext(list.files(path=args[1], pattern='matrix.csv'))
cutnames<-lapply(cutnames, function(x) gsub('matrix', '', x))
cutnames<-sort(as.character(cutnames))
files<-list.files(path=args[1], pattern='matrix.csv')
pathtofiles<-paste(args[1], files, sep='')
#clusterout<-paste(args[1], cutnames, '.cluster', sep='')
#clusterout2<-paste(args[1], 'blargh', '.cluster', sep='')
#dmats<-list()
#for (i in length(files)){
#    dtab<-read.table(files[i], header=T)
#    dmats[i]<-dist.gene(dtab)
#    }
dm_from_matrixfile <- function(filepath) {
    dist.gene(read.table(filepath, header=T, row.names=1))
}
load_table <-function(filepath){
    read.table(filepath, header=T, row.names=1)
}
dm_from_table<-function(tables){
    dist.gene(tables)
}
process<-function(pathtofiles) {
    thetable<-lapply(pathtofiles, load_table)
    distances<-lapply(thetable, dm_from_table)
    clusters <- lapply(distances, hclust)
    cuts <- lapply(clusters, function(x) cutree(x, h=0))
return(cuts)
}
cuts<-process(pathtofiles)
for(i in 1:length(cutnames)){
png(paste(args[2], cutnames[i], '.png'))
#plot(cuts[[i]], main=cutnames[i], ylab='clusternumber')
hist(cuts[[i]], main=cutnames[i])
nothing<-dev.off()
}
######
#print('check2')
#sapply((cuts), function(x) write(x, file=clusterout2))
#print('check3')
#plotted <- lapply(cuts, plot)
#print('check4')
#png(args[2])
#grid.arrange(plotted)
#n<-dev.off()
#datable<-thetable[order(sapply(thetable,length), decreasing=T)]


cutnamer<-function(cutnames, cuts){for(i in 1:length(cutnames)){
    assign(cutnames[i], cuts[[i]], envir=.GlobalEnv)}
}
cutnamer(cutnames, cuts)
strains<-names(cuts[[1]])
x<-data.table(strains)
#adds each column of the chop things and renames them to the proper name because without it they get labelled as 'i'
for(i in cutnames){
    x$i<-(get(i))
    setnames(x, 'i', i)}


permy<-as.data.table(combinations(3, 2, unlist(cutnames)))
getwallace<-function(permy){
    for(i in 1:nrow(permy))
{
#    abc<-get_abcdn(unlist(x[permy[i, V1]]), unlist(x[permy[i, V2]]))
#    walle<-wallace(unlist(abc['a']), unlist(abc['b']), unlist(abc['c']))

    argtocall1<-unlist(permy[i, 1, with=F])
    names(argtocall1)<-NULL
    argtocall2<-unlist(permy[i, 2, with=F])
    names(argtocall2)<-NULL
    adjwal<-adj_wallace(get(argtocall1), get(argtocall2))
#    write.table(walle, paste(args[1], permy[i, V1], permy[i, V2], 'wallace.csv', sep=''), row.names=F)
    write.table(adjwal, paste(args[1], permy[i, 1, with=F], permy[i, 2, with=F], 'adjwal.csv', sep=''), row.names=F)
}}
getwallace(permy)