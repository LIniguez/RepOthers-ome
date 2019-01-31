#!/bin/bash

NUMPROC=$1
FOLDERO=$2
FOLDER1=$3
FOLDER2=$4
FOLDER3=$5
GEN4ST=$6
GEN4BT=$7
MINNCOV=$8
CUT=$9

mkdir -p ${FOLDERO}

##########
# Une los archivos y convierte a SAM
########
samtools cat -o ${FOLDERO}uniq.BAM ${FOLDER1}DONE_uniq.BAM ${FOLDER2}DONE_uniq.BAM  ${FOLDER3}DONE_uniq.BAM
samtools cat -o ${FOLDERO}multiple.BAM ${FOLDER1}DONE_multiple.BAM ${FOLDER2}DONE_multiple.BAM ${FOLDER3}DONE_multiple.BAM

samtools cat -o ${FOLDERO}ALL.BAM ${FOLDERO}uniq.BAM ${FOLDERO}multiple.BAM

samtools view -@ ${NUMPROC} -b -L ${GEN4ST} ${FOLDERO}uniq.BAM > ${FOLDERO}uniq2.BAM
samtools view -@ ${NUMPROC} -b -L ${GEN4ST} ${FOLDERO}multiple.BAM > ${FOLDERO}multiple2.BAM
samtools sort -@ ${NUMPROC} -o ${FOLDERO}uniq.BAM ${FOLDERO}uniq2.BAM
samtools sort -@ ${NUMPROC} -o ${FOLDERO}multiple.BAM ${FOLDERO}multiple2.BAM

rm ${FOLDERO}uniq2.BAM ${FOLDERO}multiple2.BAM

##########
# Ve la zonas que tienen covertura para los bam de unicos y no unicos
########
bedtools genomecov -bg -ibam ${FOLDERO}uniq.BAM > ${FOLDERO}uniq_cov.bed
bedtools merge -d ${MINNCOV} -i ${FOLDERO}uniq_cov.bed > ${FOLDERO}uniq_cov_merged.bed

bedtools genomecov -bg -ibam ${FOLDERO}multiple.BAM > ${FOLDERO}multiple_cov.bed
bedtools merge -d ${MINNCOV} -i ${FOLDERO}multiple_cov.bed > ${FOLDERO}multiple_cov_merged.bed
##########
# Junta las regiones de unicos y no unicos
########
cat ${FOLDERO}multiple_cov.bed ${FOLDERO}uniq_cov.bed | sort --parallel ${NUMPROC} -V -k1,1 -k2,2n > ${FOLDERO}all_cov_sorted.bed
bedtools merge -d ${MINNCOV} -i ${FOLDERO}all_cov_sorted.bed > ${FOLDERO}all_cov_merged.bed #todas las zonas
rm ${FOLDERO}uniq_cov.bed ${FOLDERO}multiple_cov.bed  ${FOLDERO}all_cov_sorted.bed
##########
# Calcula el Coverage de las regiones
########
bedtools coverage -mean -sorted -a ${FOLDERO}uniq_cov_merged.bed -b ${FOLDERO}uniq.BAM > ${FOLDERO}uniq_4cov.bed

bedtools intersect -sorted -a ${FOLDERO}all_cov_merged.bed -b ${FOLDERO}multiple.BAM -wo >${FOLDERO}multiple_intersect_reads.bed
sort --parallel ${NUMPROC} -V -u -k 1,3 -k 7,7 ${FOLDERO}multiple_intersect_reads.bed | awk '{print $4,$5,$6}' OFS="\t" | sort --parallel ${NUMPROC} -V -k 1,1 -k2,2n > ${FOLDERO}multiple_intersect_4cov.bed
bedtools coverage -mean -sorted -a ${FOLDERO}multiple_cov_merged.bed -b ${FOLDERO}multiple_intersect_4cov.bed > ${FOLDERO}multiple_4cov.bed

bedtools intersect -wo -a ${FOLDERO}all_cov_merged.bed -b ${FOLDERO}multiple_4cov.bed ${FOLDERO}uniq_4cov.bed | awk '{ a=$7-$6;b=$8*a;print $5,$6,$7,b,a}' OFS="\t" | sort --parallel ${NUMPROC} -V -k1,1 -k2,2n > ${FOLDERO}regions_sorted_coverage.bed
bedtools merge -d ${MINNCOV} -c 4 -o sum  -i ${FOLDERO}regions_sorted_coverage.bed |awk '{a=$3-$2; b=$4/a; print $1,$2,$3,b}' OFS="\t" | awk -v CUT="$CUT" '{if (!($4 <= CUT)){ print $0;}}' > ${FOLDERO}regions_sorted_coverage_filtered.bed

rm ${FOLDERO}uniq_cov_merged.bed ${FOLDERO}multiple_cov_merged.bed ${FOLDERO}multiple_intersect_reads.bed ${FOLDERO}multiple_intersect_4cov.bed  ${FOLDERO}regions_sorted_coverage.bed
##########
# OUTPUT
########

bedtools intersect -wo -a ${FOLDERO}regions_sorted_coverage_filtered.bed -b ${FOLDERO}multiple_4cov.bed ${FOLDERO}uniq_4cov.bed | sort --parallel ${NUMPROC} -V -k 6,6 -k 7,7 > ${FOLDERO}regions_filtered_sorted.bed
rm ${FOLDERO}uniq_4cov.bed ${FOLDERO}multiple_4cov.bed

##########
# Recrea el archivo para formar la red
########
bedtools intersect -c -sorted -g ${GEN4BT} -a ${FOLDERO}regions_sorted_coverage_filtered.bed -b ${FOLDERO}uniq.BAM > ${FOLDERO}all_uniq_count.bed

bedtools intersect -g ${GEN4BT} -sorted -a ${FOLDERO}regions_sorted_coverage_filtered.bed -b ${FOLDERO}multiple.BAM -wo |sort --parallel ${NUMPROC} -V -u -k 1,3 -k 8,8 > ${FOLDERO}multiple_intersect_reads_4cov.bed

perl -e '{open(IN,"$ARGV[0]");while(<IN>){@vec=split("\t",$_);$name=$vec[0]."_".$vec[1]."_".$vec[2]; push(@{$h{$vec[7]}}, $name);}close (IN);
 foreach $k (sort keys %h){for ($cont=0;$h{$k}[$cont];$cont++){for($cont2=$cont+1;$h{$k}[$cont2];$cont2++){$out= $h{$k}[$cont]."\t".$h{$k}[$cont2];$done{$out}++;}}} undef %h;
 foreach $k (sort keys %done){print "$k\t$done{$k}\n";} }' ${FOLDERO}multiple_intersect_reads_4cov.bed > ${FOLDERO}vertex_weight_multiple.txt
perl -e '{open(RU, "$ARGV[0]");while($line=<RU>){chomp $line; @vec=split("\t",$line); $re=pop @vec;pop @vec; $name=join("_",@vec);$h{$name}=$re;}close(RU);
 open(VE, "$ARGV[1]"); while($line=<VE>){chomp $line; @vec=split("\t",$line);$w=$vec[2]+$h{$vec[0]}+$h{$vec[1]};print "$line\t$h{$vec[0]}\t$h{$vec[1]}\t$w\n";}}' ${FOLDERO}all_uniq_count.bed ${FOLDERO}vertex_weight_multiple.txt > ${FOLDERO}vertex_weight.txt
rm ${FOLDERO}uniq.BAM ${FOLDERO}vertex_weight_multiple.txt ${FOLDERO}multiple.BAM ${FOLDERO}multiple_intersect_reads_4cov.bed ${FOLDERO}all_uniq_count.bed ${FOLDERO}all_cov_merged.bed
