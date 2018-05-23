#!/bin/bash

set -e
set -u

i=0
MAX=$1
NUMPRO=$2
FOLDER=$3
FOLDERO=$4
GEN4SMT=$5
GEN4BDT=$6
BED=${FOLDER}/regions_sorted_coverage_filtered.bed


mkdir -p ${FOLDERO}

samtools sort -@ ${NUMPRO} -n -o ${FOLDERO}/${i}.BAM ${FOLDER}/ALL.BAM
samtools view -@ ${NUMPRO} ${FOLDERO}/${i}.BAM > ${FOLDERO}/${i}_temp.SAM
samtools view -H ${FOLDERO}/${i}.BAM > ${FOLDERO}/header.txt
sort --parallel ${NUMPRO} -k1,2 -V ${FOLDERO}/${i}_temp.SAM > ${FOLDERO}/${i}_sorted.SAM
cat ${FOLDERO}/header.txt ${FOLDERO}/${i}_sorted.SAM > ${FOLDERO}/${i}.SAM
rm ${FOLDERO}/${i}.BAM ${FOLDERO}/${i}_temp.SAM ${FOLDERO}/${i}_sorted.SAM


cp ${FOLDER}/vertex_weight.txt ${FOLDERO}/vertex_weight.txt
size_wei=$(echo $(wc -c ${FOLDERO}/vertex_weight.txt)| cut -f1 -d ' ')

while [ $size_wei -gt 0 ];
do
 find_communities.R ${FOLDERO}/ ${i} ${MAX} 2>/dev/null
 
 mod_txt2gtf.sh ${FOLDERO}/ ${i}
 
 ls ${FOLDERO}/*_${i}.gtf | parallel -j ${NUMPRO} -q telescope id --verbose --exp_tag ${FOLDERO}/${i}round_'{= s:_WT.+::; s:^.+${FOLDERO}/::; =}' --updated_sam ${FOLDERO}/${i}.SAM {}
 for j in $(find ${FOLDERO}/*_${i}.gtf)
 do
  k=$(echo $j| rev | cut -f 4 -d "_" | cut -f 1 -d "/"| rev)
  samtools view -@ ${NUMPRO} -b ${FOLDERO}/${i}round_${k}-updated.sam > ${FOLDERO}/${i}round_${k}-updated.bam
  samtools sort -@ ${NUMPRO} -o ${FOLDERO}/${i}round_${k}-updated_sorted.bam ${FOLDERO}/${i}round_${k}-updated.bam
  rm ${FOLDERO}/${i}round_${k}-updated.sam ${FOLDERO}/${i}round_${k}-updated.bam
 done
 
 merge_mod_bam.sh ${FOLDERO}/ ${NUMPRO} ${i} ${GEN4SMT} ${GEN4BDT} ${BED}
 i=$(($i+1))
 size_wei=$(echo $(wc -c ${FOLDERO}/vertex_weight.txt)| cut -f1 -d ' ')
done

cat ${FOLDERO}/*_0.gtf > ${FOLDERO}/ALL.gtf
cat ${FOLDERO}/header.txt ${FOLDERO}/*uniqreads.sam > ${FOLDERO}/result.sam
samtools view -@ ${NUMPRO} -b ${FOLDERO}/result.sam > ${FOLDERO}/result.bam
samtools sort -@ ${NUMPRO} -o ${FOLDERO}/result_sorted.bam ${FOLDERO}/result.bam


rm ${FOLDERO}/result.bam ${FOLDERO}/result.sam ${FOLDERO}/header.txt ${FOLDERO}/*uniqreads.sam ${FOLDERO}/*-telescope_report.tsv ${FOLDERO}/*.SAM ${FOLDERO}/*_commuity_*.gtf ${FOLDERO}/vertex_weight.txt ${FOLDERO}/*_temp.bam



