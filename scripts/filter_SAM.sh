#!/bin/sh

set -e
set -u

NUMPR=$1
FOLDER=$2
SAM=$3
SUMMARY_T=$4
EXONS=$5
RTRNACHRM=$6
MAXALL=$7

mkdir -p ${FOLDER}
#remove header
#samtools view -@ ${NUMPR} ${SAM} > ${FOLDER}/temp.sam #no es necesario solo un cp ya que te lo manda en SAM ordenado por read y best
cp ${SAM} ${FOLDER}/temp.sam
samtools view -H ${SAM} > ${FOLDER}/header.txt



>${FOLDER}/multiple.SAM
>${FOLDER}/uniq.SAM



awk -v SRR="${FOLDER}/" '{
 if(($2-256)<0){ #check if the alignment is the principal
  if($13~/XS/){  #if the alignment has XS it means it has multiple positions
   split($12,a,":");
   print $0 >> SRR"multiple.SAM"; print $1,a[3]; #print it to multiple.SAM and retain the info of the alignement score
  }else{print $0 >> SRR"uniq.SAM";}
 }else{ print $0 >> SRR"multiple.SAM";}}' OFS="\t" ${FOLDER}/temp.sam > ${FOLDER}/best_alscor.txt

#cut -f 1 ${FOLDER}/multiple.SAM | sort | uniq -c > ${FOLDER}/nummapped.txt #puede ir sin el sort ya que estan ordenados por read
cut -f 1 ${FOLDER}/multiple.SAM | uniq -c > ${FOLDER}/nummapped.txt
awk -v MAX=${MAXALL} '{ if ($1 < MAX){ print $2;}}' ${FOLDER}/nummapped.txt > ${FOLDER}/multreads_done.txt #check if the number of possible mappings is less than the maximum

perl -e '{ open(RD, "$ARGV[0]"); while($l=<RD>){chomp $l; $h{$l}=1;}
  open(BS, "$ARGV[1]"); while ($l=<BS>){ chomp $l; @vec=split("\t",$l);$bs{$vec[0]}=$vec[1];}
  open(MS, "$ARGV[2]");
  open(MSD, ">$ARGV[3]"); open(REC, ">$ARGV[4]");
  while(<MS>){		#reads SAM from multiple mapped reads
   @vec=split("\t",$_);
   if($h{$vec[0]}){	#Checks if the read is already done for mapping (it needs to have < maximum 
    @vec2=split("\:",$vec[11]);	
    if(($bs{$vec[0]}-$vec2[2])<=8){$readc{$vec[0]}++; print MSD $_;}  #checks the differences between the best alignment and the second best (only one missmatch is allowed)
   }
   elsif(($vec[1]-256)<0){ print $_;}	#if read has >= possible mappings it goes to a next round of mapping. 
  }foreach $k (keys %readc){print REC "$k\t$readc{$k}\n";} #print valid number of mappings for each read. 
  }' ${FOLDER}/multreads_done.txt ${FOLDER}/best_alscor.txt ${FOLDER}/multiple.SAM ${FOLDER}/multiple_temp.SAM ${FOLDER}/readcount.txt >  ${FOLDER}/4knext.SAM



perl -e '{ open (REC, "$ARGV[0]"); while ($l=<REC>){ chomp $l; @vec=split("\t",$l); if($vec[1]==1){$un{$vec[0]}=1;}} #A multiple mapped read can change to be uniquely mapped
 open(MS, "$ARGV[1]"); open(UNQ, ">$ARGV[2]"); 
 while (<MS>){
  @vec=split("\t",$_);
  if($un{$vec[0]}){print UNQ $_;}else{print $_;}
 }}' ${FOLDER}/readcount.txt ${FOLDER}/multiple_temp.SAM ${FOLDER}/uniq2.SAM > ${FOLDER}/multiple_done.SAM

cat ${FOLDER}/header.txt ${FOLDER}/uniq.SAM ${FOLDER}/uniq2.SAM > ${FOLDER}/unique.SAM

echo "Unique mapped reads:" >> ${SUMMARY_T}
echo $(($(wc -l ${FOLDER}/unique.SAM| cut -f 1 -d " ") - $(wc -l ${FOLDER}/header.txt| cut -f 1 -d " "))) >> ${SUMMARY_T}
echo "Multiple mapped reads: " >> ${SUMMARY_T}
echo $(cut -f 1 ${FOLDER}/multiple.SAM | sort -u | wc -l) >> ${SUMMARY_T}
echo "Reads with max mapping (k=" $MAXALL "):">> ${SUMMARY_T}
echo $(wc -l ${FOLDER}/4knext.SAM|cut -f 1 -d " ") >> ${SUMMARY_T}


cat ${FOLDER}/header.txt ${FOLDER}/multiple_done.SAM > ${FOLDER}/multiple2.SAM
mv ${FOLDER}/multiple2.SAM ${FOLDER}/multiple.SAM

cat ${FOLDER}/header.txt ${FOLDER}/4knext.SAM > ${FOLDER}/4knext2.SAM
mv ${FOLDER}/4knext2.SAM ${FOLDER}/4knext.SAM
samtools view -@ ${NUMPR} -b ${FOLDER}/4knext.SAM > ${FOLDER}/4knext.BAM
samtools fastq -n ${FOLDER}/4knext.BAM > ${FOLDER}/4knext.fastq
rm ${FOLDER}/4knext.SAM ${FOLDER}/4knext.BAM 

rm  ${FOLDER}/readcount.txt ${FOLDER}/uniq2.SAM ${FOLDER}/multiple_temp.SAM ${FOLDER}/multreads_done.txt ${FOLDER}/nummapped.txt ${FOLDER}/best_alscor.txt ${FOLDER}/temp.sam ${FOLDER}/uniq.SAM ${FOLDER}/multiple_done.SAM



samtools view -@ ${NUMPR} -h -L ${RTRNACHRM} -U ${FOLDER}/uniq_nortRNAM.SAM ${FOLDER}/unique.SAM > ${FOLDER}/rtRNAchrM.SAM
samtools view -@ ${NUMPR} -h -L ${EXONS} -U ${FOLDER}/uniq_noexons_nortRNAM.SAM ${FOLDER}/uniq_nortRNAM.SAM > ${FOLDER}/exons.SAM




echo "Reads unique mapped to exons:" >> ${SUMMARY_T}
echo $(($(wc -l ${FOLDER}/exons.SAM| cut -f 1 -d " ") - $(wc -l ${FOLDER}/header.txt| cut -f 1 -d " "))) >> ${SUMMARY_T}
echo "Reads unique mapped to rtRNA and chrM:" >> ${SUMMARY_T}
echo $(($(wc -l ${FOLDER}/rtRNAchrM.SAM| cut -f 1 -d " ") - $(wc -l ${FOLDER}/header.txt| cut -f 1 -d " "))) >> ${SUMMARY_T}
echo "Reads unique mapped left:" >> ${SUMMARY_T}
echo $(($(wc -l ${FOLDER}/uniq_noexons_nortRNAM.SAM| cut -f 1 -d " ") - $(wc -l ${FOLDER}/header.txt| cut -f 1 -d " "))) >> ${SUMMARY_T}

rm ${FOLDER}/uniq_nortRNAM.SAM ${FOLDER}/unique.SAM


samtools view -@ ${NUMPR} -L ${RTRNACHRM} -U ${FOLDER}/multiple_nortRNAM.SAM ${FOLDER}/multiple.SAM > ${FOLDER}/rtRNAchrM_multiple.SAM
cut -f 1 ${FOLDER}/rtRNAchrM_multiple.SAM | sort | uniq > ${FOLDER}/trashreads_rtRNAM.txt
cat ${FOLDER}/header.txt ${FOLDER}/multiple_nortRNAM.SAM > ${FOLDER}/multiple_nortRNAM2.SAM
mv ${FOLDER}/multiple_nortRNAM2.SAM ${FOLDER}/multiple_nortRNAM.SAM
samtools view -@ ${NUMPR} -L ${EXONS} -U ${FOLDER}/multiple_noexons_nortRNAM.SAM ${FOLDER}/multiple_nortRNAM.SAM > ${FOLDER}/exons_multiple.SAM
cut -f 1 ${FOLDER}/exons_multiple.SAM | sort | uniq > ${FOLDER}/trashreads_exons.txt


echo "Reads multiple mapped to exons:" >> ${SUMMARY_T}
echo $(wc -l ${FOLDER}/trashreads_exons.txt| cut -f 1 -d " ") >> ${SUMMARY_T}
echo "Reads multiple mapped to rtRNA and chrM:" >> ${SUMMARY_T}
echo $(wc -l ${FOLDER}/trashreads_rtRNAM.txt| cut -f 1 -d " ") >> ${SUMMARY_T}

rm ${FOLDER}/multiple_nortRNAM.SAM ${FOLDER}/exons_multiple.SAM ${FOLDER}/rtRNAchrM_multiple.SAM ${FOLDER}/multiple.SAM

cat ${FOLDER}/trashreads_exons.txt ${FOLDER}/trashreads_rtRNAM.txt >${FOLDER}/remove_multipleReads.txt
perl -e '{open(RE,"$ARGV[0]"); while($l=<RE>){chomp $l;$re{$l}=1;}close(RE);
 open(SM, "$ARGV[1]");
 while(<SM>){
  @vec=split("\t",$_);
  if(!$re{$vec[0]}){print $_;}}}' ${FOLDER}/remove_multipleReads.txt ${FOLDER}/multiple_noexons_nortRNAM.SAM >${FOLDER}/DONE_multiple_k${MAXALL}.SAM
  
rm ${FOLDER}/trashreads_exons.txt ${FOLDER}/trashreads_rtRNAM.txt ${FOLDER}/remove_multipleReads.txt ${FOLDER}/multiple_noexons_nortRNAM.SAM

echo "Multiple mapped reads left:" >>  ${SUMMARY_T}
cut -f 1 ${FOLDER}/DONE_multiple_k${MAXALL}.SAM | sort | uniq -c | awk '{print $1}' | sort | uniq -c >> ${SUMMARY_T}
cat ${FOLDER}/header.txt ${FOLDER}/DONE_multiple_k${MAXALL}.SAM > ${FOLDER}/DONE_multiple_k${MAXALL}2.SAM
mv ${FOLDER}/DONE_multiple_k${MAXALL}2.SAM ${FOLDER}/DONE_multiple_k${MAXALL}.SAM


samtools view -@ ${NUMPR} -b ${FOLDER}/exons.SAM > ${FOLDER}/exons.BAM
samtools view -@ ${NUMPR} -b ${FOLDER}/rtRNAchrM.SAM > ${FOLDER}/rtRNAchrM.BAM
samtools view -@ ${NUMPR} -b ${FOLDER}/DONE_multiple_k${MAXALL}.SAM > ${FOLDER}/DONE_multiple_k${MAXALL}.BAM
samtools view -@ ${NUMPR} -b ${FOLDER}/uniq_noexons_nortRNAM.SAM > ${FOLDER}/DONE_uniq_k${MAXALL}.BAM
samtools sort -@ ${NUMPR} -o ${FOLDER}/exons_sorted.BAM ${FOLDER}/exons.BAM
samtools sort -@ ${NUMPR} -o ${FOLDER}/rtRNAchrM_sorted.BAM ${FOLDER}/rtRNAchrM.BAM
samtools sort -@ ${NUMPR} -o ${FOLDER}/DONE_multiple_k${MAXALL}_sorted.BAM ${FOLDER}/DONE_multiple_k${MAXALL}.BAM
samtools sort -@ ${NUMPR} -o ${FOLDER}/DONE_uniq_k${MAXALL}_sorted.BAM ${FOLDER}/DONE_uniq_k${MAXALL}.BAM

rm ${FOLDER}/exons.BAM ${FOLDER}/exons.SAM ${FOLDER}/rtRNAchrM.BAM ${FOLDER}/rtRNAchrM.SAM ${FOLDER}/header.txt ${FOLDER}/uniq_noexons_nortRNAM.SAM ${FOLDER}/DONE_uniq_k${MAXALL}.BAM ${FOLDER}/DONE_multiple_k${MAXALL}.SAM ${FOLDER}/DONE_multiple_k${MAXALL}.BAM


samtools view -@ ${NUMPR} -F 256 -b -h ${FOLDER}/DONE_multiple_k${MAXALL}_sorted.BAM > ${FOLDER}/multiple.BAM
samtools fastq -n ${FOLDER}/multiple.BAM > ${FOLDER}/multiple4random.fastq
rm ${FOLDER}/multiple.BAM





