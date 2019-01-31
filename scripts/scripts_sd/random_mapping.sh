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



cat ${FLD1}multiple4random.fastq ${FLD2}multiple4random.fastq ${FLD3}multiple4random.fastq > ${FLDGR}random.fastq
bowtie2 --seed 22062018 -p ${NPRO} --very-sensitive-local --score-min L,0,1.6 --very-sensitive-local -x ${BI} -U ${FLDGR}random.fastq 2>> ${FLDG}StarDust.log | samtools view -bS -@ ${NPRO} > ${FLDGR}random.BAM
samtools cat -o ${FLDGR}random_all_rand.BAM ${FLD1}DONE_uniq.BAM ${FLD2}DONE_uniq.BAM ${FLD3}DONE_uniq.BAM ${FLDGR}random.BAM &>> ${FLDG}StarDust.log
samtools sort -@ ${NPRO} ${FLDGR}random_all_rand.BAM -o ${FLDGR}random_all_rand_sorted.BAM &>> ${FLDG}StarDust.log
rm ${FLDGR}random_all_rand.BAM ${FLDGR}random.BAM
samtools index ${FLDGR}random_all_rand_sorted.BAM
readleng=$(samtools stats -@ ${NPRO} ${FLDGR}random_all_rand_sorted.BAM | grep -P "^SN\taverage length"| awk '{print $4}')
mincov=$(samtools stats -@ ${NPRO} ${FLDGR}random_all_rand_sorted.BAM | grep -P "^RL" | sort -V -k2,2n |head -n1|cut -f 2)
numseq=$(samtools stats -@ ${NPRO} ${FLDGR}random_all_rand_sorted.BAM | grep -P "^SN\tsequences"| cut -f 3)
Rscript ${SRC}/scripts_sd/FindCoverCutoff.R ${FLDGR}random_all_rand_sorted.BAM ${readleng} 0.01 &>> ${FLDG}StarDust.log
cutoff=$(head -n 2 ${FLDGR}random_all_CoverageCutoff.txt | tail -n 1| cut -f 2) #output de FindCoverCutoff.R
rm ${FLDGR}random.fastq ${FLDGR}random*BAM* ${FLDGR}random_all_CoverageCutoff.txt
printf "Number of reads for StarDust identification:\t%d\n" $numseq>>${FLDG}summary.txt
printf "Average read length:\t%d\n" $readleng >>${FLDG}summary.txt
printf "Min read length:\t%d\n" $mincov >>${FLDG}summary.txt
printf "Coverage Cutoff:\t%f\n" $cutoff >>${FLDG}summary.txt
