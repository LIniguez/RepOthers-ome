#!/bin/bash

set -e
set -u

FLD5=$1
FLDGR=$2
FLD6=$3
NPRO=$4
CTOF=$5
GCTA=$6
GENO=$7
MNC=$8


awk '{if($5==1)print $1,$2,$3;}' OFS="\t" ${FLD5}regions_filtered_sorted.bed | sort -V -k1,1 -k2,2n | uniq > ${FLDGR}transcripts_solved_telescope.bed
perl -e '{open(IN,"$ARGV[0]"); while($l=<IN>){ chomp $l; @vec=split("\t",$l); $temp=join("_",@vec);$h{$temp}=1; }close(IN);
 open(FIL,"$ARGV[1]"); while($l=<FIL>){chomp $l; @vec=split("\t",$l); pop @vec; $temp=join("_",@vec); if(!$h{$temp}){print "$l\n";}}}' ${FLDGR}transcripts_solved_telescope.bed ${FLD5}regions_sorted_coverage_filtered.bed > ${FLDGR}transcripts_unique.bed
samtools view -b -L ${FLDGR}transcripts_unique.bed ${FLD5}ALL.BAM >${FLDGR}unique.bam
if [ $(wc -c ${FLD5}vertex_weight.txt | cut -f1 -d ' ') -gt 0 ]
then
 samtools cat -o ${FLDGR}final.bam ${FLD6}result_sorted.bam ${FLDGR}unique.bam &>> ${FLDGR}StarDust.log
 samtools sort -@ ${NPRO} -o ${FLDGR}final_sorted.bam ${FLDGR}final.bam
 bedtools genomecov -bg -ibam ${FLDGR}final_sorted.bam > ${FLDGR}regions.bed
 bedtools merge -d ${MNC} -i ${FLDGR}regions.bed > ${FLDGR}regions_merged.bed
 bedtools coverage -mean -sorted -a ${FLDGR}regions_merged.bed -b ${FLDGR}final_sorted.bam > ${FLDGR}cov.bed
 awk -v CUT="$CTOF" '{if (!($4 <= CUT)){ print $0;}}' ${FLDGR}cov.bed | sort -V -k1,1 -k2,2n > ${FLDGR}StarDust.bed
 samtools view -@ ${NPRO} -b -L ${FLDGR}StarDust.bed ${FLDGR}final_sorted.bam > ${FLDGR}StarDust.bam
 rm ${FLDGR}final_sorted.bam ${FLDGR}regions.bed ${FLDGR}regions_merged.bed ${FLDGR}cov.bed ${FLDGR}transcripts_unique.bed ${FLDGR}transcripts_solved_telescope.bed
else
 sort -V -k1,1 -k2,2n  ${FLDGR}transcripts_unique.bed > ${FLDGR}StarDust.bed
 rm ${FLDGR}transcripts_unique.bed
 mv ${FLDGR}unique.bam ${FLDGR}final.bam &>> ${FLDGR}StarDust.log
 samtools sort -@ ${NPRO} ${FLDGR}final.bam -o ${FLDGR}StarDust.bam &>> ${FLDGR}StarDust.log
fi
if [ $GCTA != "NA" ]
then
bedtools nuc -fi ${GENO} -bed ${FLDGR}StarDust.bed |awk -v THRS="$GCTA" '{if($5<THRS && $6<THRS){print $1,$2,$3,$4;} }' OFS='\t' >${FLDGR}temp.bed
mv ${FLDGR}temp.bed ${FLDGR}StarDust.bed
samtools view -@ ${NPRO} -b -L ${FLDGR}StarDust.bed ${FLDGR}StarDust.bam > ${FLDGR}temp.bam
mv ${FLDGR}temp.bam ${FLDGR}StarDust.bam
fi

awk '{a="StarDust_"NR;print $1,$2,$3,a,$4;}' OFS='\t' ${FLDGR}StarDust.bed > ${FLDGR}StarDust2.bed
mv ${FLDGR}StarDust2.bed ${FLDGR}StarDust.bed

numseq=$(samtools stats -@ ${NPRO} ${FLDGR}StarDust.bam | grep -P "^SN\tsequences"| cut -f 3)
printf "Number of Reads reported with StarDust:\t%d\n" $numseq >>${FLDGR}summary.txt
rm ${FLDGR}unique.bam ${FLDGR}final.bam
