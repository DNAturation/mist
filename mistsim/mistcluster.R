library(ape)
library(data.table)
library(tools)
library(gtools)
source('wallace-helper.R')
args<-commandArgs(trailingOnly=T)

#arg1 = tables .matrix extension
#arg2 = outfile directory

#gets the names of the saved files ready to be turned into matrices, and prepares them for future naming
cutnames<-file_path_sans_ext(list.files(path=args[1], pattern='matrix.csv'))
cutnames<-lapply(cutnames, function(x) gsub('matrix', '', x))
cutnames<-as.character(cutnames)

#gets the path of all the files
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

#function to load data into table format
load_table <-function(filepath){
    read.table(filepath, header=T, row.names=1)
}

#function to turn tables into distance matrices
dm_from_table<-function(tables){
    dist.gene(tables)
}

#function to run functions and generate clusters and cut the clusters
process<-function(pathtofiles) {
    thetable<-lapply(pathtofiles, load_table)
    distances<-lapply(thetable, dm_from_table) #order introduced onto table does not matter
    clusters <- lapply(distances, hclust)
    cuts <- lapply(clusters, function(x) cutree(x, h=0))
return(cuts)
}
cuts<-process(pathtofiles)

#creates graphs for each cluster
for(i in 1:length(cutnames)){
    png(cutnames[i])
    plot(cuts[[i]], main=cutnames[i], ylab='clusternumber')
#    hist(cuts[[i]], main=cutnames[i])
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

#sets each name of cut to call the corresponding cut
cutnamer<-function(cutnames, cuts){for(i in 1:length(cutnames)){
    assign(cutnames[i], cuts[[i]], envir=.GlobalEnv)}
}
cutnamer(cutnames, cuts)

#gets the name of all the strains. Only run on one cut due to all cuts having the same strains
strains<-names(cuts[[1]])
#adds the strain names into a data table
x<-data.table(strains)
#adds each column of the cut clusters and renames them to the proper name because without it they get labelled as 'i'
for(i in cutnames){
    x$i<-(get(i))
    setnames(x, 'i', i)}


permy<-as.data.table(combinations(length(cutnames), 2, unlist(cutnames)))
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