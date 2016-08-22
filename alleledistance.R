library(rjson)
library(dplyr)
library(argparse)
library(ggplot2)
library(reshape2)

parser <- ArgumentParser()
parser$add_argument('--jsonfile', type='character')
parser$add_argument('--outfile', type='character')
parser$add_argument('--heatmap', type='character')
args<-parser$parse_args()

jsonloader<-function(jsonfile)
  {
  json_data <- fromJSON(file=jsonfile)
  df<- bind_rows(json_data, .id = 'strain')
  nas<-sapply(df, function(x) sum(is.na(x)))
  df<-df[,c(order(nas))]
  rnas<-apply(df, 1, function(x) sum(is.na(x)))
  df<-df[c(order(rnas)),]
  return(df)
  }



thetable<-jsonloader(args$jsonfile)
write.table(thetable, args$outfile, na = '', row.names=FALSE)
thetable<-read.csv(args$outfile)
meltedtable<-melt(thetable)
ggplot(meltedtable, aes(x=strain, y=variable, fill=value)) + geom_tile()