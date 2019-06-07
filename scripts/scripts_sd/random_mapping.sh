#!/bin/bash

set -e
set -u

FLDG=$1
FLD1=$2
FLD2=$3
FLD3=$4
NPRO=$5
BI=$6
SRC=$7
SEED=$8


cat ${FLD1}multiple4random.fastq ${FLD2}multiple4random.fastq ${FLD3}multiple4random.fastq > ${FLDG}random.fastq
bowtie2 --seed 22062018 -p ${NPRO} --very-sensitive-local --score-min L,0,1.6 --very-sensitive-local -x ${BI} -U ${FLDG}random.fastq 2>> ${FLDG}StarDust.log | samtools view -bS -@ ${NPRO} > ${FLDG}random.BAM
#hisat2 --seed ${SEED} --no-unal --score-min L,0,-0.4 -p ${NPRO} --very-sensitive -x ${BI} -U ${FLDG}random.fastq 2>> ${FLDG}StarDust.log | samtools view -bS -@ ${NPRO} > ${FLDG}random.BAM
samtools view -H ${FLDG}random.BAM > ${FLDG}header.SAM
samtools cat -o ${FLDG}random_all_rand.BAM ${FLD1}DONE_uniq.BAM ${FLD2}DONE_uniq.BAM ${FLD3}DONE_uniq.BAM ${FLDG}random.BAM &>> ${FLDG}StarDust.log
samtools sort -@ ${NPRO} ${FLDG}random_all_rand.BAM -o ${FLDG}random_all_rand_sorted.BAM &>> ${FLDG}StarDust.log
rm ${FLDG}random_all_rand.BAM ${FLDG}random.BAM
samtools index ${FLDG}random_all_rand_sorted.BAM
readleng=$(samtools stats -@ ${NPRO} ${FLDG}random_all_rand_sorted.BAM | grep -P "^SN\taverage length"| awk '{print $4}')
mincov=$(samtools stats -@ ${NPRO} ${FLDG}random_all_rand_sorted.BAM | grep -P "^RL" | sort -V -k2,2n |head -n1|cut -f 2)
numseq=$(samtools stats -@ ${NPRO} ${FLDG}random_all_rand_sorted.BAM | grep -P "^SN\tsequences"| cut -f 3)
Rscript ${SRC}/scripts_sd/FindCoverCutoff.R ${FLDG}random_all_rand_sorted.BAM ${readleng} 0.01 &>> ${FLDG}StarDust.log
cutoff=$(head -n 2 ${FLDG}random_all_CoverageCutoff.txt | tail -n 1| cut -f 2) #output de FindCoverCutoff.R
rm ${FLDG}random.fastq ${FLDG}random*BAM* ${FLDG}random_all_CoverageCutoff.txt
printf "Number of reads for StarDust identification:\t%d\n" $numseq>>${FLDG}summary.txt
printf "Average read length:\t%d\n" $readleng >>${FLDG}summary.txt
printf "Min read length:\t%d\n" $mincov >>${FLDG}summary.txt
printf "Coverage Cutoff:\t%f\n" $cutoff >>${FLDG}summary.txt
