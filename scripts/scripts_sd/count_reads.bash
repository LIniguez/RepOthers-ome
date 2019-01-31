#!/bin/bash

set -e
set -u

FLDGR=$1
SRC=$2
HIS=$3
FLD4=$4
NPRO=$5
FLD1=$6
FLD2=$7
FLD3=$8
EXO=$9


awk '{b="gene_id \""$4"\"; transcript_id \""$4"\"; locus \""$4"\";"; print $1,"StarDust","transcript",$2,$3,".",".",".",b;}' OFS="\t" ${FLDGR}StarDust.bed > ${FLDGR}StarDust.gtf
Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDGR}StarDust.bam ${FLDGR}StarDust.gtf transcript gene_id ${FLDGR}StarDust FALSE &>> ${FLDGR}StarDust.log
if [ "$HIS" != "NA" ]; then
 if [ -z ${FLD4}DONE_uniq_k500.BAM ] || [ -z ${FLD4}DONE_multiple_k500.BAM ] || [ -z ${FLD4}exons_sorted.BAM ];then printf "Could not proceed spliced mapped reads were not found\n" >&2 && exit 1; fi
 samtools cat -o ${FLDGR}StarDust_splicing_temp.bam ${FLD4}DONE_uniq_k500.BAM ${FLD4}DONE_multiple_k500.BAM &>> ${FLDGR}StarDust.log
 samtools sort -@ ${NPRO} ${FLDGR}StarDust_splicing_temp.bam -o ${FLDGR}StarDust_splicing.bam &>> ${FLDGR}StarDust.log
 rm ${FLDGR}StarDust_splicing_temp.bam
 Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDGR}StarDust_splicing.bam ${FLDGR}StarDust.gtf transcript gene_id ${FLDGR}StarDust_splicing TRUE &>> ${FLDGR}StarDust.log
 samtools merge -f ${FLDGR}Exons_temp.bam ${FLD1}exons_sorted.BAM ${FLD2}exons_sorted.BAM ${FLD3}exons_sorted.BAM ${FLD4}exons_sorted.BAM &>> ${FLDGR}StarDust.log
 samtools sort -@ ${NPRO} ${FLDGR}Exons_temp.bam -o ${FLDGR}Exons.bam &>> ${FLDGR}StarDust.log
 if [ ${EXO: -4} == ".gtf" ]; then
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDGR}Exons.bam ${EXO} exon gene_id ${FLDGR}Gene TRUE &>> ${FLDGR}StarDust.log #Counting Genes and Transcript does include spliced reads
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDGR}Exons.bam ${EXO} exon transcript_id ${FLDGR}Transcript FALSE &>> ${FLDGR}StarDust.log
 fi
else
 samtools merge -f ${FLDGR}Exons_temp.bam ${FLD1}exons_sorted.BAM ${FLD2}exons_sorted.BAM ${FLD3}exons_sorted.BAM &>> ${FLDGR}StarDust.log
 samtools sort -@ ${NPRO} ${FLDGR}Exons_temp.bam -o ${FLDGR}Exons.bam &>> ${FLDGR}StarDust.log
 if [ ${EXO: -4} == ".gtf" ]; then
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDGR}Exons.bam ${EXO} exon gene_id ${FLDGR}Gene FALSE &>> ${FLDGR}StarDust.log #Counting Genes and Transcript does not include spliced reads
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDGR}Exons.bam ${EXO} exon transcript_id ${FLDGR}Transcript FALSE &>> ${FLDGR}StarDust.log
 fi
fi
 rm ${FLDGR}Exons_temp.bam
