#!/bin/bash

set -e
set -u

FLDG=$1
SRC=$2
HIS=$3
FLD4=$4
NPRO=$5
FLD1=$6
FLD2=$7
FLD3=$8
EXO=$9


awk '{b="gene_id \""$4"\"; transcript_id \""$4"\"; locus \""$4"\";"; print $1,"StarDust","transcript",$2,$3,".",".",".",b;}' OFS="\t" ${FLDG}StarDust.bed > ${FLDG}StarDust.gtf
if [ -s ${FLDG}StarDust.gtf ]; then Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDG}StarDust.bam ${FLDG}StarDust.gtf transcript gene_id ${FLDG}StarDust FALSE &>> ${FLDG}StarDust.log; fi
if [ "$HIS" != "NA" ]; then
 if [ -z ${FLD4}DONE_uniq.BAM ] || [ -z ${FLD4}DONE_multiple.BAM ] || [ -z ${FLD4}exons_sorted.BAM ];then printf "Could not proceed spliced mapped reads were not found\n" >&2 && exit 1; fi
 samtools cat -o ${FLDG}StarDust_splicing_temp.bam ${FLD4}DONE_uniq.BAM ${FLD4}DONE_multiple.BAM &>> ${FLDG}StarDust.log
 samtools sort -@ ${NPRO} ${FLDG}StarDust_splicing_temp.bam -o ${FLDG}StarDust_splicing.bam &>> ${FLDG}StarDust.log
 rm ${FLDG}StarDust_splicing_temp.bam
 if [ -s ${FLDG}StarDust.gtf ]; then Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDG}StarDust_splicing.bam ${FLDG}StarDust.gtf transcript gene_id ${FLDG}StarDust_splicing TRUE &>> ${FLDG}StarDust.log; fi;
 samtools merge -f ${FLDG}Exons_temp.bam ${FLD1}exons_sorted.BAM ${FLD2}exons_sorted.BAM ${FLD3}exons_sorted.BAM ${FLD4}exons_sorted.BAM &>> ${FLDG}StarDust.log
 samtools sort -@ ${NPRO} ${FLDG}Exons_temp.bam -o ${FLDG}Exons.bam &>> ${FLDG}StarDust.log
 if [ ${EXO: -4} == ".gtf" ]; then
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDG}Exons.bam ${EXO} exon gene_id ${FLDG}Gene TRUE &>> ${FLDG}StarDust.log #Counting Genes and Transcript does include spliced reads
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDG}Exons.bam ${EXO} exon transcript_id ${FLDG}Transcript FALSE &>> ${FLDG}StarDust.log
 fi
else
 samtools merge -f ${FLDG}Exons_temp.bam ${FLD1}exons_sorted.BAM ${FLD2}exons_sorted.BAM ${FLD3}exons_sorted.BAM &>> ${FLDG}StarDust.log
 samtools sort -@ ${NPRO} ${FLDG}Exons_temp.bam -o ${FLDG}Exons.bam &>> ${FLDG}StarDust.log
 if [ ${EXO: -4} == ".gtf" ]; then
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDG}Exons.bam ${EXO} exon gene_id ${FLDG}Gene FALSE &>> ${FLDG}StarDust.log #Counting Genes and Transcript does not include spliced reads
  Rscript ${SRC}/scripts_sd/count_transcripts.R ${FLDG}Exons.bam ${EXO} exon transcript_id ${FLDG}Transcript FALSE &>> ${FLDG}StarDust.log
 fi
fi
 rm ${FLDG}Exons_temp.bam
