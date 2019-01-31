#!/bin/bash

set -e
set -u



while getopts 'p:f:s:e:r:M:' OPTION;do
case "$OPTION" in
 p)
  NPRO=$OPTARG
 ;;
 f)
  FOLDOUT=$OPTARG
 ;;
 s)
  SMMRY=$OPTARG
 ;;
 e)
  EXO=$OPTARG
 ;;
 r)
  RT=$OPTARG
 ;;
 M)
  K=$OPTARG
 ;;
esac
done



grep -P '@SQ' ${FOLDOUT}header.txt| awk '{split($2,a,":"); split($3,b,":"); print a[2],b[2]}' OFS="\t"|sort -Vk1,1 >${FOLDOUT}header_mod.txt  #igual y se puede ordernar para que no haya problemas despues


rm ${FOLDOUT}header.txt



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
  }' ${FOLDOUT}multreads_done.txt ${FOLDOUT}best_alscor.txt ${FOLDOUT}multiple.SAM ${FOLDOUT}multiple_temp.SAM ${FOLDOUT}readcount.txt ${FOLDOUT}mult_readcount.txt >  ${FOLDOUT}4knext.SAM

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
 open(OU, ">>$ARGV[3]"); print OU "$num\n";}' ${FOLDOUT}readcount.txt ${FOLDOUT}multiple_temp.SAM ${FOLDOUT}unique.SAM ${FOLDOUT}mult_readcount.txt > ${FOLDOUT}multiple.SAM


perl -e '{ open(IN,"$ARGV[0]");$num=0;while(<IN>){@vec=split("\t",$_);
  if($_=~/^(\d+)-(\d+)\t/){$name=$1;$times=$2;
   for($i=1;$i<=$times;$i++){ $vec[0]=$name."_".$i;$o=join("\t",@vec);print $o; if(!$exists{$vec[0]}){$exists{$vec[0]}=1;$num++;}}
  }else{print $_; if(!$exists{$vec[0]} && !($vec[0] =~ /^@/)){$exists{$vec[0]}=1;$num++;} }}open(OU, ">>$ARGV[1]"); print OU "$num\n";}' ${FOLDOUT}unique.SAM ${FOLDOUT}mult_readcount.txt > ${FOLDOUT}unique_temp.SAM

mv ${FOLDOUT}unique_temp.SAM ${FOLDOUT}unique.SAM
rm ${FOLDOUT}multiple_temp.SAM


printf "Unique mapped reads:\t" >> ${SMMRY}
tail -n 1 ${FOLDOUT}mult_readcount.txt >> ${SMMRY}
printf "Multiple mapped reads:\t" >> ${SMMRY}
head -n 2 ${FOLDOUT}mult_readcount.txt | tail -n 1 >> ${SMMRY}
printf 'Reads with max mapping (k= %i):\t' $K>> ${SMMRY}
head -n 1 ${FOLDOUT}mult_readcount.txt >> ${SMMRY}


samtools view -@ ${NPRO} -b ${FOLDOUT}4knext.SAM -t ${FOLDOUT}header_mod.txt > ${FOLDOUT}4knext.BAM
samtools fastq -n ${FOLDOUT}4knext.BAM > ${FOLDOUT}4knext.fastq 2>/dev/null
rm ${FOLDOUT}4knext.SAM ${FOLDOUT}4knext.BAM ${FOLDOUT}readcount.txt ${FOLDOUT}multreads_done.txt ${FOLDOUT}best_alscor.txt ${FOLDOUT}mult_readcount.txt

if [ $RT != "NA.bed" ]
then
 samtools view -@ ${NPRO} -b -h -L ${RT} -U ${FOLDOUT}uniq_nortRNAM.BAM -t ${FOLDOUT}header_mod.txt ${FOLDOUT}unique.SAM > ${FOLDOUT}rtRNAchrM.BAM
 samtools view -@ ${NPRO} -b -h -L ${EXO} -U ${FOLDOUT}uniq_noexons_nortRNAM.BAM ${FOLDOUT}uniq_nortRNAM.BAM > ${FOLDOUT}exons_sorted.BAM
 printf "Reads unique mapped to rtRNA and chrM:" >> ${SMMRY}
 echo $(samtools stats ${FOLDOUT}rtRNAchrM.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SMMRY}
 rm ${FOLDOUT}uniq_nortRNAM.BAM
 samtools view -@ ${NPRO} -L ${RT} -U ${FOLDOUT}multiple_nortRNAM.SAM -t ${FOLDOUT}header_mod.txt ${FOLDOUT}multiple.SAM |cut -f 1 | sort --parallel ${NPRO} -u > ${FOLDOUT}trashreads.txt
 samtools view -@ ${NPRO} -L ${EXO} -U ${FOLDOUT}multiple_noexons_nortRNAM.SAM -t ${FOLDOUT}header_mod.txt ${FOLDOUT}multiple_nortRNAM.SAM | cut -f 1 | sort --parallel ${NPRO} -u >> ${FOLDOUT}trashreads.txt
 rm ${FOLDOUT}multiple_nortRNAM.SAM ${FOLDOUT}rtRNAchrM.BAM
else
 samtools view -@ ${NPRO} -b -h -L ${EXO} -U ${FOLDOUT}uniq_noexons_nortRNAM.BAM -t ${FOLDOUT}header_mod.txt ${FOLDOUT}unique.SAM > ${FOLDOUT}exons_sorted.BAM
 samtools view -@ ${NPRO} -L ${EXO} -U ${FOLDOUT}multiple_noexons_nortRNAM.SAM -t ${FOLDOUT}header_mod.txt ${FOLDOUT}multiple.SAM | cut -f 1 | sort --parallel ${NPRO} -u > ${FOLDOUT}trashreads.txt
fi

printf "Reads unique mapped to exons:\t" >> ${SMMRY}
echo $(samtools stats ${FOLDOUT}exons_sorted.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SMMRY}
printf "Reads multiple mapped eliminated:\t" >> ${SMMRY}
echo $(wc -l ${FOLDOUT}trashreads.txt| cut -f 1 -d " ") >> ${SMMRY}
printf "Reads unique mapped left:\t" >> ${SMMRY}
echo $(samtools stats ${FOLDOUT}uniq_noexons_nortRNAM.BAM | grep -P '^SN\traw total' | cut -f 3) >> ${SMMRY}
rm ${FOLDOUT}multiple.SAM ${FOLDOUT}unique.SAM

perl -e '{open(RE1,"$ARGV[0]"); while($l=<RE1>){chomp $l;$re{$l}=1;}close(RE1);
 open(SM, "$ARGV[1]"); open (CNT, ">$ARGV[2]");
 while(<SM>){
  @vec=split("\t",$_);
  if(!$re{$vec[0]}){print $_; $h{$vec[0]}++;}
 }foreach $read (keys %h){print CNT "$read\t$h{$read}\n";}
 }' ${FOLDOUT}trashreads.txt ${FOLDOUT}multiple_noexons_nortRNAM.SAM ${FOLDOUT}reads_passed.txt >${FOLDOUT}DONE_multiple.SAM

rm ${FOLDOUT}trashreads.txt ${FOLDOUT}multiple_noexons_nortRNAM.SAM

printf "Multiple mapped reads left:\n" >>  ${SMMRY}
cut -f 2 ${FOLDOUT}reads_passed.txt | sort -n --parallel ${NPRO} | uniq -c >> ${SMMRY}

samtools view -@ ${NPRO} -b ${FOLDOUT}DONE_multiple.SAM -t ${FOLDOUT}header_mod.txt > ${FOLDOUT}DONE_multiple.BAM
mv ${FOLDOUT}uniq_noexons_nortRNAM.BAM ${FOLDOUT}DONE_uniq.BAM

rm ${FOLDOUT}header_mod.txt  ${FOLDOUT}DONE_multiple.SAM  ${FOLDOUT}reads_passed.txt

samtools view -@ ${NPRO} -F 256 -b -h ${FOLDOUT}DONE_multiple.BAM > ${FOLDOUT}multiple.BAM
samtools fastq -n ${FOLDOUT}multiple.BAM > ${FOLDOUT}multiple4random.fastq 2>/dev/null
rm ${FOLDOUT}multiple.BAM
