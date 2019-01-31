#!/bin/bash

set -e
set -u

I=0
MAX=$1
WTM=$2
NPRO=$3
FLD5=$4
FOLDOUT=$5
GEN4ST=$6
GEN4BT=$7
SRC=$8
SEED=$9

samtools view -h -@ ${NPRO} ${FLD5}ALL.BAM > ${FOLDOUT}${I}.SAM
samtools view -H ${FOLDOUT}${I}.SAM >${FOLDOUT}header.txt


cp ${FLD5}vertex_weight.txt ${FOLDOUT}vertex_weight.txt
size_wei=$(echo $(wc -c ${FOLDOUT}vertex_weight.txt)| cut -f1 -d ' ')

while [ $size_wei -gt 0 ];
do
 Rscript ${SRC}/find_communities.R ${FOLDOUT} ${I} ${MAX} ${WTM} ${SEED} 2>/dev/null

 sh ${SRC}/mod_txt2gtf.sh ${FOLDOUT} ${I}
 ls ${FOLDOUT}*_${I}.gtf | parallel -j ${NPRO} -q telescope assign --outdir ${FOLDOUT} --tempdir ${FOLDOUT} --quiet --exp_tag ${I}round_'{= s:_WT.+::; s:^.+${FOLDOUT}/::; =}' --updated_sam ${FOLDOUT}${I}.SAM {}
 for j in $(find ${FOLDOUT}*_${I}.gtf)
 do
  k=$(echo $j| rev | cut -f 4 -d "_" | cut -f 1 -d "/"| rev)
  samtools sort -@ ${NPRO} -o ${FOLDOUT}${I}round_${k}-updated_sorted.bam ${FOLDOUT}${I}round_${k}-updated.bam
  rm ${FOLDOUT}${I}round_${k}-updated.bam ${FOLDOUT}${I}round_${k}-checkpoint.npz ${FOLDOUT}${I}round_${k}-tmp_tele.bam ${FOLDOUT}${I}round_${k}-other.bam
 done

 sh ${SRC}/merge_mod_bam.sh ${FOLDOUT} ${NPRO} ${I} ${GEN4ST} ${GEN4BT} ${FLD5}regions_sorted_coverage_filtered.bed
 I=$(($I+1))
 size_wei=$(echo $(wc -c ${FOLDOUT}vertex_weight.txt)| cut -f1 -d ' ')
done

cat ${FOLDOUT}*_0.gtf > ${FOLDOUT}ALL.gtf
cat ${FOLDOUT}header.txt ${FOLDOUT}*uniqreads.sam > ${FOLDOUT}result.sam
samtools view -@ ${NPRO} -b ${FOLDOUT}result.sam > ${FOLDOUT}result.bam
samtools sort -@ ${NPRO} -o ${FOLDOUT}result_sorted.bam ${FOLDOUT}result.bam


rm ${FOLDOUT}result.bam ${FOLDOUT}result.sam ${FOLDOUT}header.txt ${FOLDOUT}*uniqreads.sam ${FOLDOUT}*-telescope_report.tsv ${FOLDOUT}*.SAM ${FOLDOUT}*_commuity_*.gtf ${FOLDOUT}vertex_weight.txt ${FOLDOUT}*_telescope_res*
