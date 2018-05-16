#!/usr/bin/env Rscript


args<-commandArgs(trailingOnly=T)
bamfile<-args[1]
annotation<-args[2]
feature<-args[3]
type<-args[4]
outcou<-paste0(args[5],"_count.txt")
outspl<-paste0(args[5],"_junctions.txt")
splicing<-args[6]

library(Rsubread)

fc<-featureCounts(bamfile, annot.ext=annotation, isGTFAnnotationFile=T, GTF.featureType=feature, GTF.attrType=type, allowMultiOverlap=T, juncCounts=T, countMultiMappingReads=T)
write.table(fc$counts, file=outcou, sep="\t",quote=F)

if(as.logical(splicing)){
junc<-fc$counts_junction[(!is.na((fc$counts_junction)$PrimaryGene)),]
write.table(junc, file=outspl, sep="\t",quote=F)
}