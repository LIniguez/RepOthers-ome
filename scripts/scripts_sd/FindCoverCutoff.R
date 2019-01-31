#!/usr/bin/env Rscript


library(transcriptR)
args <- commandArgs(trailingOnly=TRUE)

check <- strsplit(args[1], "_")[[1]]
samp <- paste(check[1:(length(check)-2)],collapse="_")
samp_o<-paste0(name,"_CoverageCutoff.txt")
tds <- constructTDS(file = args[1], fragment.size = as.numeric(args[2]),unique = FALSE)
estimateBackground(tds, fdr.cutoff = as.numeric(args[3]))
write.table(tds@coverageCutoff,file=samp_o, sep="\t",quote=F)
