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


check_sam(){
 >${1}unique.SAM; >${1}best_alscor.txt; >${1}header.txt
 awk -v FOLD="$1" '{
 if($1~/^@/) {print $0 >> FOLD"header.txt";} #Remove header sequencs
 if(($2-256)<0){ #check if the alignment is the principal
  if($13~/XS/){  #if the alignment has XS it means it has multiple positions in Bowtie2
   split($12,a,":");print $0 ; print $1,a[3] >> FOLD"best_alscor.txt"; #print it to multiple.SAM and retain the info of the alignement score
  }else if($13~/ZS/){  #if the alignment has XS it means it has multiple positions in Bowtie2
   split($12,a,":"); print $0 ; print $1,a[3] >> FOLD"best_alscor.txt"; #print it to multiple.SAM and retain the info of the alignement score
  }else{print $0 >> FOLD"unique.SAM";}}else{ print $0 ;}}' OFS="\t"
}
check_mult(){
 FOLD=$1; MAX=$2
 perl -pe '{@vec=split("\t",$_);$h{$vec[0]}++;END{open(OUT, ">$OUT");
 foreach $read(keys %h){if( $h{$read} < $MAX ){ print OUT "$read\n";}}}}' -s -- -OUT=${FOLD}multreads_done.txt -MAX=${MAX}
}

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
