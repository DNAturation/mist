library(ape)
library(data.table)
library(tools)
args<-commandArgs(trailingOnly=T)

#1 = core_calls.csv
#2 = chopsym directory
#3 = outdirectory + number for file
df<-read.csv(args[1])
DT <- data.table(df)
#colnames(DT)[names(DT)=='X'] <- 'Strain'
passnames<-file_path_sans_ext(list.files(path=args[2], pattern='.fasta'))
DT2<-DT[, passnames, with=F]
DT2<-cbind(DT2, DT[, .(Strain=X)])
#DT2<-DT[, .(Strain=X), by=X]

#moves the strain names to the first column
setcolorder(DT2, c('Strain', passnames))

#genedists <- dist.gene(DT2)
write.table(DT2, args[3], row.names=F)