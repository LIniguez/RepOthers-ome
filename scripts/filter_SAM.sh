#!/bin/sh

set -e
set -u

# NUMPR=$1
# FOLDER=$2
# SAM=$3
# SUMMARY_T=$4
# EXONS=$5
# RTRNACHRM=$6
# MAXALL=$7

NUMPR=$1
FOLDER=$2
SUMMARY_T=$3
EXONS=$4
RTRNACHRM=$5
MAXALL=$6


#mkdir -p ${FOLDER}
#remove header
#samtools view -@ ${NUMPR} ${SAM} > ${FOLDER}/temp.sam #no es necesario solo un cp ya que te lo manda en SAM ordenado por read y best
#mv ${SAM} ${FOLDER}/temp.sam
#samtools view -H ${FOLDER}/temp.sam > ${FOLDER}/header.txt
grep -P '@SQ' ${FOLDER}/header.txt| awk '{split($2,a,":"); split($3,b,":"); print a[2],b[2]}' OFS="\t" >${FOLDER}/header_mod.txt
rm ${FOLDER}/header.txt


# >${FOLDER}/multiple.SAM
# >${FOLDER}/uniq.SAM
# awk -v SRR="${FOLDER}/" '{
 # if($1~/^@/) {next;} #Remove header sequencs
 # if(($2-256)<0){ #check if the alignment is the principal
  # if($13~/XS/){  #if the alignment has XS it means it has multiple positions
   # split($12,a,":");
   # print $0 >> SRR"multiple.SAM"; print $1,a[3]; #print it to multiple.SAM and retain the info of the alignement score
  # }else{print $0 >> SRR"uniq.SAM";}
 # }else{ print $0 >> SRR"multiple.SAM";}}' OFS="\t" ${FOLDER}/temp.sam > ${FOLDER}/best_alscor.txt
# parallel --pipepart --tmpdir ${FOLDER}/ --block -1 -j ${NUMPR} -a ${FOLDER}/temp.sam -q awk -v SRR="${FOLDER}/{%}_" '{
 # if($1~/^@/) {next;} #Remove header sequencs
 # if(($2-256)<0){ #check if the alignment is the principal
  # if($13~/XS/){  #if the alignment has XS it means it has multiple positions
   # split($12,a,":");
   # print $0 >> SRR"multiple.SAM"; print $1,a[3]; #print it to multiple.SAM and retain the info of the alignement score
  # }else{print $0 >> SRR"uniq.SAM";}
 # }else{ print $0 >> SRR"multiple.SAM";}}' OFS="\t"  > ${FOLDER}/best_alscor.txt
# if ls ${FOLDER}/*_uniq.SAM >/dev/null 2>&1 
# then
 # cat ${FOLDER}/*_uniq.SAM > ${FOLDER}/unique.SAM
 # rm ${FOLDER}/*_uniq.SAM
# else
 # > ${FOLDER}/unique.SAM
# fi
# if ls ${FOLDER}/*_multiple.SAM >/dev/null 2>&1 
# then
 # cat ${FOLDER}/*_multiple.SAM > ${FOLDER}/multiple.SAM
 # rm ${FOLDER}/*_multiple.SAM
# else
 # > ${FOLDER}/multiple.SAM
# fi
#cut -f 1 ${FOLDER}/multiple.SAM | sort --parallel ${NUMPR} | uniq -c > ${FOLDER}/nummapped.txt
#parallel --pipepart --block -1 -j ${NUMPR} -a ${FOLDER}/nummapped.txt -q awk -v MAX=${MAXALL} '{ if ($1 < MAX){ print $2;}}'  > ${FOLDER}/multreads_done.txt
# perl -e '{open(IN,"$ARGV[0]");
 # while(<IN>){
 # @vec=split("\t",$_);$h{$vec[0]}++;}
 # foreach $read(keys %h){
 # if( $h{$read} < $ARGV[1] ){ print "$read\n";}}}' ${FOLDER}/multiple.SAM ${MAXALL}> ${FOLDER}/multreads_done.txt
#cut -f 1 ${FOLDER}/multiple.SAM | uniq -c > ${FOLDER}/nummapped.txt
#awk -v MAX=${MAXALL} '{ if ($1 < MAX){ print $2;}}' ${FOLDER}/nummapped.txt > ${FOLDER}/multreads_done.txt #check if the number of possible mappings is less than the maximum

perl -e '{ open(RD, "$ARGV[0]"); while($l=<RD>){chomp $l; $h{$l}=1;} #reads with less number of mappings than maximum
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
 open(MS, "$ARGV[1]"); open(UNQ, ">>$ARGV[2]"); 
 while (<MS>){
  @vec=split("\t",$_);
  if($un{$vec[0]}){print UNQ $_;}else{print $_;}
 }}' ${FOLDER}/readcount.txt ${FOLDER}/multiple_temp.SAM ${FOLDER}/unique.SAM > ${FOLDER}/multiple.SAM


perl -e '{ open(IN,"$ARGV[0]");while(<IN>){if($_=~/^(\d+)-(\d+)\t/){$name=$1;$times=$2; @vec=split("\t",$_);for($i=1;$i<=$times;$i++){ $vec[0]=$name."_".$i;$o=join("\t",@vec);print $o;}}else{print $_;}}}' ${FOLDER}/unique.SAM > ${FOLDER}/unique_temp.SAM &
perl -e '{ open(IN,"$ARGV[0]");while(<IN>){if($_=~/^(\d+)-(\d+)\t/){$name=$1;$times=$2; @vec=split("\t",$_);for($i=1;$i<=$times;$i++){ $vec[0]=$name."_".$i;$o=join("\t",@vec);print $o;}}else{print $_;}}}' ${FOLDER}/multiple.SAM > ${FOLDER}/multiple_temp.SAM &
wait
mv ${FOLDER}/unique_temp.SAM ${FOLDER}/unique.SAM
mv ${FOLDER}/multiple_temp.SAM ${FOLDER}/multiple.SAM

echo "Unique mapped reads:" >> ${SUMMARY_T}
echo $(wc -l ${FOLDER}/unique.SAM| cut -f 1 -d " ") >> ${SUMMARY_T}
echo "Multiple mapped reads: " >> ${SUMMARY_T}
echo $(cut -f 1 ${FOLDER}/multiple.SAM | sort --parallel ${NUMPR} -u | wc -l) >> ${SUMMARY_T}
echo "Reads with max mapping (k=" $MAXALL "):">> ${SUMMARY_T}
echo $(wc -l ${FOLDER}/4knext.SAM|cut -f 1 -d " ") >> ${SUMMARY_T}

samtools view -@ ${NUMPR} -b ${FOLDER}/4knext.SAM -t ${FOLDER}/header_mod.txt > ${FOLDER}/4knext.BAM
samtools fastq -n ${FOLDER}/4knext.BAM > ${FOLDER}/4knext.fastq 2>/dev/null
rm ${FOLDER}/4knext.SAM ${FOLDER}/4knext.BAM ${FOLDER}/readcount.txt ${FOLDER}/multreads_done.txt ${FOLDER}/best_alscor.txt

sed -e '/^\s*$/d' ${FOLDER}/unique.SAM > ${FOLDER}/uniq2.SAM
mv ${FOLDER}/uniq2.SAM ${FOLDER}/unique.SAM
samtools view -@ ${NUMPR} -b -h -L ${RTRNACHRM} -U ${FOLDER}/uniq_nortRNAM.BAM -t ${FOLDER}/header_mod.txt ${FOLDER}/unique.SAM > ${FOLDER}/rtRNAchrM.BAM
samtools view -@ ${NUMPR} -b -h -L ${EXONS} -U ${FOLDER}/uniq_noexons_nortRNAM.BAM ${FOLDER}/uniq_nortRNAM.BAM > ${FOLDER}/exons.BAM

echo "Reads unique mapped to exons:" >> ${SUMMARY_T}
echo $(samtools stats ${FOLDER}/exons.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SUMMARY_T}
echo "Reads unique mapped to rtRNA and chrM:" >> ${SUMMARY_T}
echo $(samtools stats ${FOLDER}/rtRNAchrM.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SUMMARY_T}
echo "Reads unique mapped left:" >> ${SUMMARY_T}
echo $(samtools stats ${FOLDER}/uniq_noexons_nortRNAM.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SUMMARY_T}

rm ${FOLDER}/uniq_nortRNAM.BAM ${FOLDER}/unique.SAM

samtools view -@ ${NUMPR} -L ${RTRNACHRM} -U ${FOLDER}/multiple_nortRNAM.SAM -t ${FOLDER}/header_mod.txt ${FOLDER}/multiple.SAM > ${FOLDER}/rtRNAchrM_multiple.SAM
cut -f 1 ${FOLDER}/rtRNAchrM_multiple.SAM | sort --parallel ${NUMPR} -u > ${FOLDER}/trashreads_rtRNAM.txt
samtools view -@ ${NUMPR} -L ${EXONS} -U ${FOLDER}/multiple_noexons_nortRNAM.SAM -t ${FOLDER}/header_mod.txt ${FOLDER}/multiple_nortRNAM.SAM > ${FOLDER}/exons_multiple.SAM
cut -f 1 ${FOLDER}/exons_multiple.SAM | sort --parallel ${NUMPR} -u > ${FOLDER}/trashreads_exons.txt
	
echo "Reads multiple mapped to exons:" >> ${SUMMARY_T}
echo $(wc -l ${FOLDER}/trashreads_exons.txt| cut -f 1 -d " ") >> ${SUMMARY_T}
echo "Reads multiple mapped to rtRNA and chrM:" >> ${SUMMARY_T}
echo $(wc -l ${FOLDER}/trashreads_rtRNAM.txt| cut -f 1 -d " ") >> ${SUMMARY_T}

rm ${FOLDER}/multiple_nortRNAM.SAM ${FOLDER}/exons_multiple.SAM ${FOLDER}/rtRNAchrM_multiple.SAM ${FOLDER}/multiple.SAM

perl -e '{open(RE1,"$ARGV[0]"); while($l=<RE1>){chomp $l;$re{$l}=1;}close(RE1);
 open(RE2,"$ARGV[1]"); while($l=<RE2>){chomp $l;$re{$l}=1;}close(RE2);
 open(SM, "$ARGV[2]"); open (CNT, ">$ARGV[3]");
 while(<SM>){
  @vec=split("\t",$_);
  if(!$re{$vec[0]}){print $_; $h{$vec[0]}++;}
 }foreach $read (keys %h){print CNT "$read\t$h{$read}\n";}
 }' ${FOLDER}/trashreads_exons.txt ${FOLDER}/trashreads_rtRNAM.txt ${FOLDER}/multiple_noexons_nortRNAM.SAM ${FOLDER}/reads_passed.txt >${FOLDER}/DONE_multiple_k${MAXALL}.SAM
  
rm ${FOLDER}/trashreads_exons.txt ${FOLDER}/trashreads_rtRNAM.txt ${FOLDER}/multiple_noexons_nortRNAM.SAM 

echo "Multiple mapped reads left:" >>  ${SUMMARY_T}
cut -f 2 ${FOLDER}/reads_passed.txt | sort --parallel ${NUMPR} | uniq -c >> ${SUMMARY_T}

samtools view -@ ${NUMPR} -b ${FOLDER}/DONE_multiple_k${MAXALL}.SAM -t ${FOLDER}/header_mod.txt > ${FOLDER}/DONE_multiple_k${MAXALL}.BAM
samtools sort -@ ${NUMPR} -o ${FOLDER}/exons_sorted.BAM ${FOLDER}/exons.BAM 2>/dev/null
samtools sort -@ ${NUMPR} -o ${FOLDER}/rtRNAchrM_sorted.BAM ${FOLDER}/rtRNAchrM.BAM 2>/dev/null
samtools sort -@ ${NUMPR} -o ${FOLDER}/DONE_multiple_k${MAXALL}_sorted.BAM ${FOLDER}/DONE_multiple_k${MAXALL}.BAM 2>/dev/null
samtools sort -@ ${NUMPR} -o ${FOLDER}/DONE_uniq_k${MAXALL}_sorted.BAM ${FOLDER}/uniq_noexons_nortRNAM.BAM 2>/dev/null

rm ${FOLDER}/exons.BAM ${FOLDER}/rtRNAchrM.BAM ${FOLDER}/header_mod.txt ${FOLDER}/uniq_noexons_nortRNAM.BAM ${FOLDER}/DONE_multiple_k${MAXALL}.SAM ${FOLDER}/DONE_multiple_k${MAXALL}.BAM ${FOLDER}/reads_passed.txt

samtools view -@ ${NUMPR} -F 256 -b -h ${FOLDER}/DONE_multiple_k${MAXALL}_sorted.BAM > ${FOLDER}/multiple.BAM
samtools fastq -n ${FOLDER}/multiple.BAM > ${FOLDER}/multiple4random.fastq 2>/dev/null
rm ${FOLDER}/multiple.BAM





