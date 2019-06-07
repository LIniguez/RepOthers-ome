#!/usr/bin/env Rscript

library('derfinder')

#threshold<-round(getTotalMapped(rawFile="rail-rna_out/coverage_bigwigs/simulation.bw")/1000000) * 0.25
threshold<-7
datos<-fullCoverage(files="rail-rna_out/coverage_bigwigs/simulation.bw", chrs=paste0('chr',c(1:21,"X","Y")),  totalMapped = 1, targetSize = 1)
regionMat<-regionMatrix(datos, cutoff = threshold, L = 100)
for (i in paste0('chr',c(1:21,"X","Y"))){rtracklayer::export.bed(regionMat[[i]]$regions, paste0(i,"_Railderfinder.bed"))}
for (i in paste0('chr',c(1:21,"X","Y"))){write.table(round(regionMat[[i]]$coverageMatrix),paste0(i,"_Railderfinder_count.txt"),sep="\t",row.names = F, col.names = F)}
system("cat *_Railderfinder.bed > Railderfinder.bed")
system("rm *_Railderfinder.bed")
system("cat *_Railderfinder_count.txt > Railderfinder_count.txt")
system("rm *_Railderfinder_count.txt")
system("paste Railderfinder.bed Railderfinder_count.txt >Output_Railderfinder_summarized.bed")