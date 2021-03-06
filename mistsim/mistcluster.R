#may need to add in something that compares arg startchop to how many files there are, as if the input is greater or
# equal to the number of files, it currently crashes, so exact numbers are required (if greater, cannot find files. if
#equal, calling dist.gene on the last set with nothing to compare to raises an error)

library(ape)
library(data.table)
library(tools)
library(gtools)
library(parallel)
library(compiler)
library(argparse)
source('wallace-helper.R')
#args<-commandArgs(trailingOnly=T)

#arg1 = core_calls.csv
#arg2 = chopsym directory
#arg3 = 'chopped.csv'
#arg4 = outpath, 'misty/sim/
#arg5 = outcuts, 'misty/cuts/
#arg6 = startchops
parser <- ArgumentParser()
parser$add_argument('corecalls', type='character')
parser$add_argument('--chopsyms', type='character')
parser$add_argument('--chopped', type='character')
parser$add_argument('--outpath', type='character')
parser$add_argument('--outcuts', type='character')
parser$add_argument('--startchops', nargs='+')
args<-parser$parse_args()


dg<-cmpfun(dist.gene)

dm_from_matrixfile <- function(filepath) {
    dg(read.table(filepath, header=T, row.names=1))
}

#function to load data into table format
load_table <-function(filepath){
    read.table(filepath, header=T, row.names=1)
}

#function to turn tables into distance matrices
#nNOTE: this takes a long time to run, try to multicore it?
dm_from_table<-function(tables){
    dg(tables)
}

#grabs the corecalls file and stores it into memory
corecalls<-data.table(read.csv(paste(args$corecalls)))
#gets gene names in order of worst to best
genes<-as.list(data.table(read.csv(paste(args$outpath, unlist(args$startchops)[1], args$chopped, sep=''), header=F))[, 1, with=F])
names(genes)[1]<-unlist(args$startchops)[1]
#gets set names of genes assigned to the cut they belong to
for (i in 2:length(unlist(args$startchops))){ ###can this be vectored instead of looped? Opening a csv file, difficult
    temp <- as.list(data.table(read.csv(paste(args$outpath, args$startchops[i], args$chopped, sep=''), header=F))[, 1, with=F])
    genes<-c(genes, temp)
    names(genes)[i]<-unlist(args$startchops[i])
}


#if (length(unlist(args$startchops)) == 2){###can this be vectored instead of looped?
#    temp <- as.list(data.table(read.csv(paste(args$outpath, args$startchops[2], args$chopped, sep='')))[, 1, with=F])
#    genes<-c(genes, temp)
#    names(genes)[2]<-unlist(args$startchops[2])
#    }
#
#if (length(unlist(args$startchops)) == 1){
#    next}

#assume get gene[i] = list of a bunch of genes
clustermaker<-function(genelists, nameofgenes){
    glist<-unlist(genelists)
    newDT<-corecalls[, glist, with=F]
    distances<-dm_from_table(newDT)
    clusters<-hclust(distances)
    assign(paste('cuts', nameofgenes, sep=''), cutree(clusters, h=0))
    write.table(get(paste('cuts', nameofgenes, sep='')), paste(args$outcuts, 'cut', nameofgenes, '.csv', sep=''))
    return(get(paste('cuts', nameofgenes, sep='')))
}
cutsmatrix<-mcmapply(clustermaker, genes, names(genes), mc.cores=60)



#gets the name of all the strains.
strains<-corecalls[, X]
#adds the strain names into a data table
wallacetable<-data.table(strains)
#adds each column of the cut clusters and renames them to the proper name because without it they get labelled as 'i'
for(i in names(genes)){
    wallacetable$i<-cutsmatrix[, toString(i)]
    setnames(wallacetable, 'i', paste('cuts', i, sep=''))
    }

#get adjusted wallace values for each chop compared against every other chop, DEFAULT VERSION
#permy<-as.data.table(combinations(length(cutnames), 2, unlist(cutnames)))
#getwallace<-function(permy){
#    for(i in 1:nrow(permy))
#{
##    abc<-get_abcdn(unlist(x[permy[i, V1]]), unlist(x[permy[i, V2]]))
##    walle<-wallace(unlist(abc['a']), unlist(abc['b']), unlist(abc['c']))
#
#    argtocall1<-unlist(permy[i, 1, with=F])
#    names(argtocall1)<-NULL
#    argtocall2<-unlist(permy[i, 2, with=F])
#    names(argtocall2)<-NULL
#    adjwal<-adj_wallace(get(argtocall1), get(argtocall2))
##    write.table(walle, paste(args[1], permy[i, V1], permy[i, V2], 'wallace.csv', sep=''), row.names=F)
#    write.table(adjwal, paste(args[1], permy[i, 1, with=F], permy[i, 2, with=F], 'adjwal.csv', sep=''), row.names=F)
#}}
#getwallace(permy)




#compares all to a single reference for adjwallace,
# also writes out the adjusted wallace results to a csv file (can be removed)
ordered<-sort(as.integer(names(genes)))
temp<-adj_wallace(unlist(wallacetable[, paste('cuts', ordered[1], sep=''), with=F]),
    unlist(wallacetable[, paste('cuts', ordered[1], sep=''), with=F]))
wallacelist<-temp['Adjusted_Wallace_B_vs_A']
names(wallacelist)<-paste('cuts', ordered[1], sep='')
for(i in 2:length(ordered)){

    argtocall1<-unlist(wallacetable[, paste('cuts', ordered[1], sep=''), with=F])
    names(argtocall1)<-NULL
    argtocall2<-unlist(wallacetable[, paste('cuts', ordered[i], sep=''), with=F])
    names(argtocall2)<-NULL
    adjwal<-adj_wallace(argtocall1, argtocall2)
    write.table(adjwal, paste(args$outcuts, ordered[1], '_', ordered[i], 'adjwal.csv', sep=''), row.names=F)
    temp<-adjwal['Adjusted_Wallace_B_vs_A']
    names(temp) = paste('cuts', ordered[i], sep='')
    wallacelist<-append(wallacelist, temp)

}

wallaceplot<-cbind(as.integer(ordered), wallacelist)
png('wallaceplot.png')
plot(wallaceplot)
nothing<-dev.off()
