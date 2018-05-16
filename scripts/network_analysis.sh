#!/bin/bash

set -e
set -u

i=0
MAX=$1
FOLDER=$2
FOLDERO=$3
GEN4SMT=$4
BED=${FOLDER}/regions_sorted_coverage_filtered.bed


mkdir -p ${FOLDERO}

samtools view -@ 8 -h ${FOLDER}/ALL.BAM > ${FOLDERO}/${i}.SAM
cp ${FOLDER}/vertex_weight.txt ${FOLDERO}/vertex_weight.txt
size_wei=$(echo $(wc -c ${FOLDERO}/vertex_weight.txt)| cut -f1 -d ' ')

while [ $size_wei -gt 0 ];
do
 find_communities.R ${FOLDERO}/ ${i} ${MAX}
 
 mod_txt2gtf.sh ${FOLDERO}/ ${i}
 
 for j in $(find ${FOLDERO}/*_${i}.gtf)
 do
  k=$(echo $j| rev | cut -f 2 -d "_" | rev)
  telescope id --verbose --exp_tag ${FOLDERO}/${i}round_${k} --updated_sam ${FOLDERO}/${i}.SAM ${j}
  samtools view -@ 8 -b ${FOLDERO}/${i}round_${k}-updated.sam > ${FOLDERO}/${i}round_${k}-updated.bam
  samtools sort -o ${FOLDERO}/${i}round_${k}-updated_sorted.bam ${FOLDERO}/${i}round_${k}-updated.bam
  rm ${FOLDERO}/${i}round_${k}-updated.sam ${FOLDERO}/${i}round_${k}-updated.bam
 done
 
 merge_mod_bam.sh ${FOLDERO}/ ${i} ${GEN4SMT} ${BED}
 i=$(($i+1))
 size_wei=$(echo $(wc -c ${FOLDERO}/vertex_weight.txt)| cut -f1 -d ' ')

done

cat ${FOLDERO}/*_0.gtf > ${FOLDERO}/ALL.gtf
cat ${FOLDERO}/header.txt ${FOLDERO}/*uniqreads.sam > ${FOLDERO}/result.sam
samtools view -b ${FOLDERO}/result.sam > ${FOLDERO}/result.bam
samtools sort -o ${FOLDERO}/result_sorted.bam ${FOLDERO}/result.bam


rm ${FOLDERO}/result.bam ${FOLDERO}/result.sam ${FOLDERO}/header.txt ${FOLDERO}/*uniqreads.sam ${FOLDERO}/*-telescope_report.tsv ${FOLDERO}/*.SAM ${FOLDERO}/*_commuity_*.gtf ${FOLDERO}/vertex_weight.txt ${FOLDERO}/*_temp.bam



