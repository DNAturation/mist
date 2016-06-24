library(ape)
library(data.table)
library(tools)
args<-commandArgs(trailingOnly=T)

#1 = core_calls.csv
#2 = chopsym directory
#3 = outdirectory + number for file
df<-read.csv(args[1])
DT <- data.table(df)
names<-file_path_sans_ext(list.files(path=args[2], pattern='.fasta'))
DT2<-DT[, .(Strain=X), by=names]
setcolorder(DT2, c('Strain', names))
#genedists <- dist.gene(DT2)
write.table(DT2, args[3], row.names=F)