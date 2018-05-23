awk '{a=$1"_"$2"_"$3;b="gene_id \""a"\";"; print $1,"repeatsome","gene",$2,$3,".",".",".",b;}' OFS="\t" GENES_Prostate.bed >genes_repeatsome.gtf
bedtools intersect -a GENES_Prostate.bed -b filtered_isoforms.bed -wao | awk '{a=$1"_"$2"_"$3;c=a"."$4"_"$5"_"$6;b="gene_id \""a"\"; transcript_id \""c"\";"; print $4,"repeatsome","transcript",$5,$6,".",".",".",b;}' OFS="\t" > transcripts_repeatsome.gtf

cat genes.gtf transcripts_repeatsome.gtf | sort -V -k1,1 -k4,4n >REPEATSOME.gtf
