samp=Lassa_project_Human

folds=/home/ubuntu/efs/lassa/RepOthers/*
length=101

cat ${folds}/Rep*bed | sort -V -k1,1 -k2,2n > Allpossible_transcripts_RepOthes.bed

bedtools merge -i Allpossible_transcripts_RepOthes.bed > ${samp}_RepOther_genes.bed

bedtools intersect -a ${samp}_RepOther_genes.bed -b Allpossible_transcripts_RepOthes.bed -wao > Intersection_genes_isoforms.bed

find_isoforms.pl Intersection_genes_isoforms.bed ${length} > filtered_isoforms_${samp}.bed

awk '{a=$1"_"$2"_"$3;b="gene_id \""a"\";"; print $1,"repeatsome","gene",$2,$3,".",".",".",b;}' OFS="\t" ${samp}_RepOther_genes.bed >${samp}_RepOther_genes.gtf

bedtools intersect -a ${samp}_RepOther_genes.bed -b filtered_isoforms_${samp}.bed -wao | awk '{a=$1"_"$2"_"$3;c=a"."$4"_"$5"_"$6;b="gene_id \""a"\"; transcript_id \""c"\";"; print $4,"repeatsome","transcript",$5,$6,".",".",".",b;}' OFS="\t" > transcripts_RepOthers-ome_${samp}.gtf

cat ${samp}_RepOther_genes.gtf transcripts_RepOthers-ome_${samp}.gtf | sort -V -k1,1 -k4,4n >${samp}_RepOthers.gtf

rm Allpossible_transcripts_RepOthes.bed Intersection_genes_isoforms.bed ${samp}_RepOther_genes.bed transcripts_RepOthers-ome_${samp}.gtf

htseq-count -f bam -s no -a 0 -q -t transcript ${folds}/RepOthers.bam ${samp}_RepOthers.gtf >htseq-count-RepOthers_${samp}.txt


#############
#
#R script
#
##############
samp<-"prostate". #just the name of the output

read_multiple<-function(x, name, col){
    lapply(lfile, function(x) {
    z <-read.delim(x,header=T,skip=1)[,col]
    colnames(z)<-c(name,strsplit(x,split="/")[[1]][1])    #You need to check this because in my case I have my outputs saved in a folder named after the sample. Example :SRR497707/XXXYYY- -telescope_report.tsv but you need to check this
    z
})}.   #this function is for reading multiple files

lfile<-list.files(pattern=' -telescope_report.tsv', recursive=T) #find all files with that name
datatemp <- read_multiple(lfile, name="transcript", col=c(1,4)) # Read the file and extract the column 3 is for final_count and 4 the final_conf
datos<-Reduce(function(x,y) {merge(x,y, by="transcript", all=T)}, datatemp) #merge all tables into a single one
datos[is.na(datos)]<-0
file1<-paste0("Telescope_results_table_",samp,".txt")
write.table(datos, file=file1, quote=F, sep='\t', row.names=F)
