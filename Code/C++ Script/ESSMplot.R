#1个输入，ESS matrix的路径
args<-commandArgs(T)
ESSM<-read.csv(args[1], header=TRUE, row.names=1, sep="\t")

library(ggplot2)
library(pheatmap)
library(RColorBrewer)

breaksList = seq(0,1,by = 0.01)
Heatmapplot <- pheatmap(ESSM, cluster_rows = FALSE, cluster_cols = FALSE, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)

outputpath <- paste(args[1],".Heatmap.png",sep = "")

ggsave(outputpath,Heatmapplot,width = 3, height = 8)
