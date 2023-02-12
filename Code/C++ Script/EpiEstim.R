#1个输入，DailyCases的路径
args<-commandArgs(T)
Daily<-read.csv(args[1], header=TRUE, row.names=1, sep="\t")

library(devtools)
library(EpiEstim)

Rdf <- data.frame(matrix(ncol = 0,nrow = nrow(Daily)-7))

for(i in 1:ncol(Daily))
{
res <- estimate_R(Daily[,i], method = "parametric_si",config = make_config(list(mean_si = 5.5, std_si = 4.5)))
Rdf <- cbind(Rdf,res$R$'Mean(R)')
}

names(Rdf) <- names(Daily)

outputpath <- paste(args[1],".RE",sep = "")
write.table(Rdf,outputpath,sep = "\t")