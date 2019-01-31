#!/bin/bash

NPRO=$1
FOLDOUT=$2
FLD1=$3
FLD2=$4
FLD3=$5
GEN4ST=$6
GEN4BT=$7
MNC=$8
CUT=$9


##########
# Une los archivos y convierte a SAM
########
samtools cat -o ${FOLDOUT}uniq.BAM ${FLD1}DONE_uniq.BAM ${FLD2}DONE_uniq.BAM  ${FLD3}DONE_uniq.BAM
samtools cat -o ${FOLDOUT}multiple.BAM ${FLD1}DONE_multiple.BAM ${FLD2}DONE_multiple.BAM ${FLD3}DONE_multiple.BAM

samtools cat -o ${FOLDOUT}ALL.BAM ${FOLDOUT}uniq.BAM ${FOLDOUT}multiple.BAM

samtools view -@ ${NPRO} -b -L ${GEN4ST} ${FOLDOUT}uniq.BAM > ${FOLDOUT}uniq2.BAM
samtools view -@ ${NPRO} -b -L ${GEN4ST} ${FOLDOUT}multiple.BAM > ${FOLDOUT}multiple2.BAM
samtools sort -@ ${NPRO} -o ${FOLDOUT}uniq.BAM ${FOLDOUT}uniq2.BAM
samtools sort -@ ${NPRO} -o ${FOLDOUT}multiple.BAM ${FOLDOUT}multiple2.BAM

rm ${FOLDOUT}uniq2.BAM ${FOLDOUT}multiple2.BAM

##########
# Ve la zonas que tienen covertura para los bam de unicos y no unicos
########
bedtools genomecov -bg -ibam ${FOLDOUT}uniq.BAM > ${FOLDOUT}uniq_cov.bed
bedtools merge -d ${MNC} -i ${FOLDOUT}uniq_cov.bed > ${FOLDOUT}uniq_cov_merged.bed

bedtools genomecov -bg -ibam ${FOLDOUT}multiple.BAM > ${FOLDOUT}multiple_cov.bed
bedtools merge -d ${MNC} -i ${FOLDOUT}multiple_cov.bed > ${FOLDOUT}multiple_cov_merged.bed
##########
# Junta las regiones de unicos y no unicos
########
cat ${FOLDOUT}multiple_cov.bed ${FOLDOUT}uniq_cov.bed | sort --parallel ${NPRO} -V -k1,1 -k2,2n > ${FOLDOUT}all_cov_sorted.bed
bedtools merge -d ${MNC} -i ${FOLDOUT}all_cov_sorted.bed > ${FOLDOUT}all_cov_merged.bed #todas las zonas
rm ${FOLDOUT}uniq_cov.bed ${FOLDOUT}multiple_cov.bed  ${FOLDOUT}all_cov_sorted.bed
##########
# Calcula el Coverage de las regiones
########
bedtools coverage -mean -sorted -a ${FOLDOUT}uniq_cov_merged.bed -b ${FOLDOUT}uniq.BAM > ${FOLDOUT}uniq_4cov.bed

bedtools intersect -sorted -a ${FOLDOUT}all_cov_merged.bed -b ${FOLDOUT}multiple.BAM -wo >${FOLDOUT}multiple_intersect_reads.bed
sort --parallel ${NPRO} -V -u -k 1,3 -k 7,7 ${FOLDOUT}multiple_intersect_reads.bed | awk '{print $4,$5,$6}' OFS="\t" | sort --parallel ${NPRO} -V -k 1,1 -k2,2n > ${FOLDOUT}multiple_intersect_4cov.bed
bedtools coverage -mean -sorted -a ${FOLDOUT}multiple_cov_merged.bed -b ${FOLDOUT}multiple_intersect_4cov.bed > ${FOLDOUT}multiple_4cov.bed

bedtools intersect -wo -a ${FOLDOUT}all_cov_merged.bed -b ${FOLDOUT}multiple_4cov.bed ${FOLDOUT}uniq_4cov.bed | awk '{ a=$7-$6;b=$8*a;print $5,$6,$7,b,a}' OFS="\t" | sort --parallel ${NPRO} -V -k1,1 -k2,2n > ${FOLDOUT}regions_sorted_coverage.bed
bedtools merge -d ${MNC} -c 4 -o sum  -i ${FOLDOUT}regions_sorted_coverage.bed |awk '{a=$3-$2; b=$4/a; print $1,$2,$3,b}' OFS="\t" | awk -v CUT="$CUT" '{if (!($4 <= CUT)){ print $0;}}' > ${FOLDOUT}regions_sorted_coverage_filtered.bed

rm ${FOLDOUT}uniq_cov_merged.bed ${FOLDOUT}multiple_cov_merged.bed ${FOLDOUT}multiple_intersect_reads.bed ${FOLDOUT}multiple_intersect_4cov.bed  ${FOLDOUT}regions_sorted_coverage.bed
##########
# OUTPUT
########

bedtools intersect -wo -a ${FOLDOUT}regions_sorted_coverage_filtered.bed -b ${FOLDOUT}multiple_4cov.bed ${FOLDOUT}uniq_4cov.bed | sort --parallel ${NPRO} -V -k 6,6 -k 7,7 > ${FOLDOUT}regions_filtered_sorted.bed
rm ${FOLDOUT}uniq_4cov.bed ${FOLDOUT}multiple_4cov.bed

##########
# Recrea el archivo para formar la red
########
bedtools intersect -c -sorted -g ${GEN4BT} -a ${FOLDOUT}regions_sorted_coverage_filtered.bed -b ${FOLDOUT}uniq.BAM > ${FOLDOUT}all_uniq_count.bed

bedtools intersect -g ${GEN4BT} -sorted -a ${FOLDOUT}regions_sorted_coverage_filtered.bed -b ${FOLDOUT}multiple.BAM -wo |sort --parallel ${NPRO} -V -u -k 1,3 -k 8,8 > ${FOLDOUT}multiple_intersect_reads_4cov.bed

perl -e '{open(IN,"$ARGV[0]");while(<IN>){@vec=split("\t",$_);$name=$vec[0]."_".$vec[1]."_".$vec[2]; push(@{$h{$vec[7]}}, $name);}close (IN);
 foreach $k (sort keys %h){for ($cont=0;$h{$k}[$cont];$cont++){for($cont2=$cont+1;$h{$k}[$cont2];$cont2++){$out= $h{$k}[$cont]."\t".$h{$k}[$cont2];$done{$out}++;}}} undef %h;
 foreach $k (sort keys %done){print "$k\t$done{$k}\n";} }' ${FOLDOUT}multiple_intersect_reads_4cov.bed > ${FOLDOUT}vertex_weight_multiple.txt
perl -e '{open(RU, "$ARGV[0]");while($line=<RU>){chomp $line; @vec=split("\t",$line); $re=pop @vec;pop @vec; $name=join("_",@vec);$h{$name}=$re;}close(RU);
 open(VE, "$ARGV[1]"); while($line=<VE>){chomp $line; @vec=split("\t",$line);$w=$vec[2]+$h{$vec[0]}+$h{$vec[1]};print "$line\t$h{$vec[0]}\t$h{$vec[1]}\t$w\n";}}' ${FOLDOUT}all_uniq_count.bed ${FOLDOUT}vertex_weight_multiple.txt > ${FOLDOUT}vertex_weight.txt
rm ${FOLDOUT}uniq.BAM ${FOLDOUT}vertex_weight_multiple.txt ${FOLDOUT}multiple.BAM ${FOLDOUT}multiple_intersect_reads_4cov.bed ${FOLDOUT}all_uniq_count.bed ${FOLDOUT}all_cov_merged.bed
