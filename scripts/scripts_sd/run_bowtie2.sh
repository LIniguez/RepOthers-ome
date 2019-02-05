#!/bin/bash

set -e
set -u

SEED=$1
UNAL=$2
NPRO=$3
K=$4
BI=$5
FASQ=$6
FOLDOUT=$7
FLDG=$8
RT=$9
SRC=${10}


. ${SRC}/scripts_sd/functions.sh

if [ ${UNAL} != "NA" ]; then
  bowtie2 --seed ${SEED} --un-gz ${UNAL} --no-unal --score-min L,0,1.6 -p ${NPRO} -k ${K} --very-sensitive-local -x ${BI} -U ${FASQ} 2>> ${FLDG}StarDust.log | check_sam ${FOLDOUT} | check_mult ${FOLDOUT} ${K} > ${FOLDOUT}multiple.SAM
else
  if [ ${FASQ: -5} == "fasta"  ]
  then
    bowtie2 --seed ${SEED} -f --no-unal --score-min L,0,1.6 -p ${NPRO} -k ${K} --very-sensitive-local -x ${BI} -U ${FASQ} 2>> ${FLDG}StarDust.log | check_sam ${FOLDOUT} | check_mult ${FOLDOUT} ${K} > ${FOLDOUT}multiple.SAM
  else
    bowtie2 --seed ${SEED} --no-unal --score-min L,0,1.6 -p ${NPRO} -k ${K} --very-sensitive-local -x ${BI} -U ${FASQ} 2>> ${FLDG}StarDust.log | check_sam ${FOLDOUT} | check_mult ${FOLDOUT} ${K} > ${FOLDOUT}multiple.SAM
  fi
fi

bash ${SRC}/scripts_sd/filter_SAM.sh -p ${NPRO} -f ${FOLDOUT} -s ${FLDG}summary.txt -e ${FLDG}exons_anotation.bed -r ${RT} -M ${K} &>> ${FLDG}StarDust.log
