#!/bin/bash

set -e
set -u



while getopts 'p:f:s:e:r:M:' OPTION;do
case "$OPTION" in
 p)
  NUMPR=$OPTARG
 ;;
 f)
  FOLDER=$OPTARG
 ;;
 s)
  SUMMARY_T=$OPTARG
 ;;
 e)
  EXONS=$OPTARG
 ;;
 r)
  RTRNACHRM=$OPTARG
 ;;
 M)
  MAXALL=$OPTARG
 ;;
esac
done



grep -P '@SQ' ${FOLDER}header.txt| awk '{split($2,a,":"); split($3,b,":"); print a[2],b[2]}' OFS="\t"|sort -Vk1,1 >${FOLDER}header_mod.txt  #igual y se puede ordernar para que no haya problemas despues


rm ${FOLDER}header.txt



perl -e '{ open(RD, "$ARGV[0]"); while($l=<RD>){chomp $l; $h{$l}=1;} close(RD); #reads with less number of mappings than maximum
  open(BS, "$ARGV[1]"); while ($l=<BS>){ chomp $l; @vec=split("\t",$l);$bs{$vec[0]}=$vec[1];} close(BS);
  open(MS, "$ARGV[2]");
  open(MSD, ">$ARGV[3]");
  open(REC, ">$ARGV[4]");$num=0;
  while(<MS>){		#reads SAM from multiple mapped reads
   @vec=split("\t",$_);
   if($h{$vec[0]}){	#Checks if the read is already done for mapping (it needs to have < maximum )
    @vec2=split("\:",$vec[11]);
    if(($bs{$vec[0]}-$vec2[2])<=8){$readc{$vec[0]}++; print MSD $_;}  #checks the differences between the best alignment and the second best (only one missmatch is allowed)
   }
   elsif(($vec[1]-256)<0){ print $_;$num++;}	#if read has >= possible mappings it goes to a next round of mapping.
  }foreach $k (keys %readc){print REC "$k\t$readc{$k}\n"; } open (OU, ">$ARGV[5]");print OU "$num\n"; #print valid number of mappings for each read.
  }' ${FOLDER}multreads_done.txt ${FOLDER}best_alscor.txt ${FOLDER}multiple.SAM ${FOLDER}multiple_temp.SAM ${FOLDER}readcount.txt ${FOLDER}mult_readcount.txt >  ${FOLDER}4knext.SAM

perl -e '{ open (REC, "$ARGV[0]"); while ($l=<REC>){ chomp $l; @vec=split("\t",$l); if($vec[1]==1){$un{$vec[0]}=1;}} #A multiple mapped read can change to be uniquely mapped
 open(MS, "$ARGV[1]"); open(UNQ, ">>$ARGV[2]"); $num=0;
 while (<MS>){
  @vec=split("\t",$_);
  if($un{$vec[0]}){print UNQ $_;}
  else{
   if($_=~/^(\d+)-(\d+)\t/){$name=$1;$times=$2; @vec=split("\t",$_);
    for($i=1;$i<=$times;$i++){ $vec[0]=$name."_".$i;$o=join("\t",@vec);print $o; if(!$exists{$vec[0]}){$exists{$vec[0]}=1;$num++;}}
   }else{print $_; if(!$exists{$vec[0]}){$exists{$vec[0]}=1;$num++;} }} #changes for reads collapsed into a single read (after -k 4)
 }
 open(OU, ">>$ARGV[3]"); print OU "$num\n";}' ${FOLDER}readcount.txt ${FOLDER}multiple_temp.SAM ${FOLDER}unique.SAM ${FOLDER}mult_readcount.txt > ${FOLDER}multiple.SAM


perl -e '{ open(IN,"$ARGV[0]");$num=0;while(<IN>){@vec=split("\t",$_);
  if($_=~/^(\d+)-(\d+)\t/){$name=$1;$times=$2;
   for($i=1;$i<=$times;$i++){ $vec[0]=$name."_".$i;$o=join("\t",@vec);print $o; if(!$exists{$vec[0]}){$exists{$vec[0]}=1;$num++;}}
  }else{print $_; if(!$exists{$vec[0]} && !($vec[0] =~ /^@/)){$exists{$vec[0]}=1;$num++;} }}open(OU, ">>$ARGV[1]"); print OU "$num\n";}' ${FOLDER}unique.SAM ${FOLDER}mult_readcount.txt > ${FOLDER}unique_temp.SAM

mv ${FOLDER}unique_temp.SAM ${FOLDER}unique.SAM
rm ${FOLDER}multiple_temp.SAM


echo "Unique mapped reads:" >> ${SUMMARY_T}
tail -n 1 ${FOLDER}mult_readcount.txt >> ${SUMMARY_T}
echo "Multiple mapped reads: " >> ${SUMMARY_T}
head -n 2 ${FOLDER}mult_readcount.txt | tail -n 1 >> ${SUMMARY_T}
echo "Reads with max mapping (k=" $MAXALL "):">> ${SUMMARY_T}
head -n 1 ${FOLDER}mult_readcount.txt >> ${SUMMARY_T}


samtools view -@ ${NUMPR} -b ${FOLDER}4knext.SAM -t ${FOLDER}header_mod.txt > ${FOLDER}4knext.BAM
samtools fastq -n ${FOLDER}4knext.BAM > ${FOLDER}4knext.fastq 2>/dev/null
rm ${FOLDER}4knext.SAM ${FOLDER}4knext.BAM ${FOLDER}readcount.txt ${FOLDER}multreads_done.txt ${FOLDER}best_alscor.txt ${FOLDER}mult_readcount.txt

if [ $RTRNACHRM != "NA.bed" ]
then
 samtools view -@ ${NUMPR} -b -h -L ${RTRNACHRM} -U ${FOLDER}uniq_nortRNAM.BAM -t ${FOLDER}header_mod.txt ${FOLDER}unique.SAM > ${FOLDER}rtRNAchrM.BAM
 samtools view -@ ${NUMPR} -b -h -L ${EXONS} -U ${FOLDER}uniq_noexons_nortRNAM.BAM ${FOLDER}uniq_nortRNAM.BAM > ${FOLDER}exons_sorted.BAM
 echo "Reads unique mapped to rtRNA and chrM:" >> ${SUMMARY_T}
 echo $(samtools stats ${FOLDER}rtRNAchrM.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SUMMARY_T}
 rm ${FOLDER}uniq_nortRNAM.BAM
 samtools view -@ ${NUMPR} -L ${RTRNACHRM} -U ${FOLDER}multiple_nortRNAM.SAM -t ${FOLDER}header_mod.txt ${FOLDER}multiple.SAM |cut -f 1 | sort --parallel ${NUMPR} -u > ${FOLDER}trashreads.txt
 samtools view -@ ${NUMPR} -L ${EXONS} -U ${FOLDER}multiple_noexons_nortRNAM.SAM -t ${FOLDER}header_mod.txt ${FOLDER}multiple_nortRNAM.SAM | cut -f 1 | sort --parallel ${NUMPR} -u >> ${FOLDER}trashreads.txt
 rm ${FOLDER}multiple_nortRNAM.SAM ${FOLDER}rtRNAchrM.BAM
else
 samtools view -@ ${NUMPR} -b -h -L ${EXONS} -U ${FOLDER}uniq_noexons_nortRNAM.BAM -t ${FOLDER}header_mod.txt ${FOLDER}unique.SAM > ${FOLDER}exons_sorted.BAM
 samtools view -@ ${NUMPR} -L ${EXONS} -U ${FOLDER}multiple_noexons_nortRNAM.SAM -t ${FOLDER}header_mod.txt ${FOLDER}multiple.SAM | cut -f 1 | sort --parallel ${NUMPR} -u > ${FOLDER}trashreads.txt
fi

echo "Reads unique mapped to exons:" >> ${SUMMARY_T}
echo $(samtools stats ${FOLDER}exons_sorted.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SUMMARY_T}
echo "Reads multiple mapped eliminated:" >> ${SUMMARY_T}
echo $(wc -l ${FOLDER}trashreads.txt| cut -f 1 -d " ") >> ${SUMMARY_T}
echo "Reads unique mapped left:" >> ${SUMMARY_T}
echo $(samtools stats ${FOLDER}uniq_noexons_nortRNAM.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SUMMARY_T}
rm ${FOLDER}multiple.SAM ${FOLDER}unique.SAM

perl -e '{open(RE1,"$ARGV[0]"); while($l=<RE1>){chomp $l;$re{$l}=1;}close(RE1);
 open(SM, "$ARGV[1]"); open (CNT, ">$ARGV[2]");
 while(<SM>){
  @vec=split("\t",$_);
  if(!$re{$vec[0]}){print $_; $h{$vec[0]}++;}
 }foreach $read (keys %h){print CNT "$read\t$h{$read}\n";}
 }' ${FOLDER}trashreads.txt ${FOLDER}multiple_noexons_nortRNAM.SAM ${FOLDER}reads_passed.txt >${FOLDER}DONE_multiple_k${MAXALL}.SAM

rm ${FOLDER}trashreads.txt ${FOLDER}multiple_noexons_nortRNAM.SAM

echo "Multiple mapped reads left:" >>  ${SUMMARY_T}
cut -f 2 ${FOLDER}reads_passed.txt | sort --parallel ${NUMPR} | uniq -c >> ${SUMMARY_T}

samtools view -@ ${NUMPR} -b ${FOLDER}DONE_multiple_k${MAXALL}.SAM -t ${FOLDER}header_mod.txt > ${FOLDER}DONE_multiple_k${MAXALL}.BAM
mv ${FOLDER}uniq_noexons_nortRNAM.BAM ${FOLDER}DONE_uniq_k${MAXALL}.BAM

rm ${FOLDER}header_mod.txt  ${FOLDER}DONE_multiple_k${MAXALL}.SAM  ${FOLDER}reads_passed.txt

samtools view -@ ${NUMPR} -F 256 -b -h ${FOLDER}DONE_multiple_k${MAXALL}.BAM > ${FOLDER}multiple.BAM
samtools fastq -n ${FOLDER}multiple.BAM > ${FOLDER}multiple4random.fastq 2>/dev/null
rm ${FOLDER}multiple.BAM
