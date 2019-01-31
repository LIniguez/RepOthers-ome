#!/bin/bash

set -e
set -u

FLDGR=$1
RANDFASQ=${FLDGR}random.fastq
OUT_RAND=${FLDGR}random
FLD1=$2
FLD2=$3
FLD3=$4
NPRO=$5
BI=$6
SRC=$7




cat ${FLD1}multiple4random.fastq ${FLD2}multiple4random.fastq ${FLD3}multiple4random.fastq > ${RANDFASQ}
bowtie2 --seed 22062018 -p ${NPRO} --very-sensitive-local --score-min L,0,1.6 --very-sensitive-local -x ${BI} -U ${RANDFASQ} 2>> ${FLDGR}StarDust.log | samtools view -bS -@ ${NPRO} > ${OUT_RAND}.BAM
samtools cat -o ${OUT_RAND}_all_rand.BAM ${FLD1}DONE_uniq.BAM ${FLD2}DONE_uniq.BAM ${FLD3}DONE_uniq.BAM ${OUT_RAND}.BAM &>> ${FLDGR}StarDust.log
samtools sort -@ ${NPRO} ${OUT_RAND}_all_rand.BAM -o ${OUT_RAND}_all_rand_sorted.BAM &>> ${FLDGR}StarDust.log
rm ${OUT_RAND}_all_rand.BAM ${OUT_RAND}.BAM
samtools index ${OUT_RAND}_all_rand_sorted.BAM
readleng=$(samtools stats -@ ${NPRO} ${OUT_RAND}_all_rand_sorted.BAM | grep -P "^SN\taverage length"| awk '{print $4}')
mincov=$(samtools stats -@ ${NPRO} ${OUT_RAND}_all_rand_sorted.BAM | grep -P "^RL" | sort -V -k2,2n |head -n1|cut -f 2)
numseq=$(samtools stats -@ ${NPRO} ${OUT_RAND}_all_rand_sorted.BAM | grep -P "^SN\tsequences"| cut -f 3)
Rscript ${SRC}/scripts_sd/FindCoverCutoff.R ${OUT_RAND}_all_rand_sorted.BAM ${readleng} 0.01 &>> ${FLDGR}StarDust.log
cutoff=$(head -n 2 ${OUT_RAND}_all_CoverageCutoff.txt | tail -n 1| cut -f 2) #output de FindCoverCutoff.R
rm ${RANDFASQ} ${OUT_RAND}*BAM* ${OUT_RAND}_all_CoverageCutoff.txt
printf "Number of reads for StarDust identification:\t%d\n" $numseq>>${FLDGR}summary.txt
printf "Average read length:\t%d\n" $readleng >>${FLDGR}summary.txt
printf "Min read length:\t%d\n" $mincov >>${FLDGR}summary.txt
printf "Coverage Cutoff:\t%f\n" $cutoff >>${FLDGR}summary.txt
