#!/bin/bash

set -e
set -u
numpro=1
startin=4
folder_gral=./
bow_index=""
rtRNA_MChr=""
exons=""
hi_index=""
fastq1=""
fastq2=""
fastq=""

while getopts 'p:b:r:e:o:i:1:2:U:hS:' OPTION;do
 case "$OPTION" in
 p)
  numpro=$OPTARG
  ;;
 b)
  bow_index=$OPTARG
  ;;
 r)
  rtRNA_MChr=$OPTARG
  ;;
 e)
  exons=$OPTARG
  ;;
 o)
  folder_gral=$OPTARG
  ;;
 i)
  hi_index=$OPTARG
  ;;
 h)
  echo "Aqui viene una ayuda"
  exit 1
  ;;
 1)
  fastq1=$OPTARG
  ;;
 2)
  fastq2=$OPTARG
  ;;
 U)
  fastq=$OPTARG
  ;;
 S)
  startin=$OPTARG
  ;;
 ?)
  echo "script usage: $(basename $0) [-l] [-h] [-a somevalue]" >&2 #mensaje de ayuda para correr comando en caso de que haya un parametro que no venga al caso
  exit 1
  ;;
 esac
done

if [ -z "$bow_index" ] || [ -z "$rtRNA_MChr" ] || [ -z "$exons" ] || [ -z "$hi_index" ] 
then
 echo "Not all annotations/DB needed" && exit 1
fi

if [ -z "$fastq1" ] || [ -z "$fastq2" ]
then
 if [ -z "$fastq" ]
 then
  echo "Missing fastq file(s)" && exit 1
 else
  paired=FALSE	
 fi
else
 paired=TRUE
 fastq=${folder_gral}/test_joined.fastq
fi




check_sam(){
 >${1}/unique.SAM
 >${1}/best_alscor.txt
 >${1}/header.txt
 awk -v FOLD="$1" '{
 if($1~/^@/) {print $0 >> FOLD"/header.txt";} #Remove header sequencs
 if(($2-256)<0){ #check if the alignment is the principal
  if($13~/XS/){  #if the alignment has XS it means it has multiple positions
   split($12,a,":");
   print $0 ; print $1,a[3] >> FOLD"/best_alscor.txt"; #print it to multiple.SAM and retain the info of the alignement score
  }else{print $0 >> FOLD"/unique.SAM";}
 }else{ print $0 ;}}' OFS="\t"  
}
check_mult(){
 FOLD=$1
 MAX=$2
 perl -pe '{
  @vec=split("\t",$_);$h{$vec[0]}++;
  END{open(OUT, ">$OUT");
   foreach $read(keys %h){
   if( $h{$read} < $MAX ){ print OUT "$read\n";}}}}' -s -- -OUT=${FOLD}/multreads_done.txt -MAX=${MAX}
}
unali=${folder_gral}/unalign_temp.gz
unali2=${folder_gral}/unalign.bz2
folder1=${folder_gral}/mapping_4
folder2=${folder_gral}/mapping_100
folder3=${folder_gral}/mapping_500
folder4=${folder_gral}/mapping_splicesites
folder5=${folder_gral}/transcripts
folder6=${folder_gral}/network
randfastq=${folder_gral}/for_random.fq
out_rand=${folder_gral}/random


if [ $startin == "4" ] 
then 
 mkdir -p ${folder_gral}/
 if [ $paired == "TRUE" ];
 then
  perl -e '{open(IN,$ARGV[0]); while(<IN>){@vec=split(" ",$_);$l=$vec[0]."_".$ARGV[1]."\n"; $a=<IN>;$b=<IN>;$c=<IN>;print "$l$a+\n$c";}}' ${fastq1} 1 >${fastq}
  perl -e '{open(IN,$ARGV[0]); while(<IN>){@vec=split(" ",$_);$l=$vec[0]."_".$ARGV[1]."\n"; $a=<IN>;$b=<IN>;$c=<IN>;print "$l$a+\n$c";}}' ${fastq2} 2 >>${fastq}
 fi
 awk '{if($3 == "exon")print $1,$4,$5;}' OFS="\t" ${exons} | sort --parallel ${numpro} -V -k1,1 -k2,2n |bedtools merge -i stdin >${folder_gral}/exons_anotation.bed
 >${folder_gral}/summary.txt
 >${folder_gral}/RepOthers-ome.log
 mkdir -p ${folder1}/
 echo "Mapping"
 echo " Bowtie2 -k4"
 #First Bowtie2 with a small number of possible alignments
 bowtie2 --seed 22062018 --un-gz ${unali} --no-unal --score-min L,0,1.6 -p ${numpro} -k 4 --very-sensitive-local -x ${bow_index} -U ${fastq} 2>> ${folder_gral}/RepOthers-ome.log | check_sam ${folder1} | check_mult ${folder1} 4 > ${folder1}/multiple.SAM
 filter_SAM.sh ${numpro} ${folder1} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 4 &>> ${folder_gral}/RepOthers-ome.log
 fastx_collapser -i ${folder1}/4knext.fastq -o ${folder1}/4knext.fasta
 rm ${folder1}/4knext.fastq

fi


if [ $startin == "4" ] || [ $startin == "100" ] 
then 
 mkdir -p ${folder2}/
 echo " Bowtie2 -k100"
 #Second round
 bowtie2 --seed 22062018 --no-unal --score-min L,0,1.6 -f -p ${numpro} -k 100 --very-sensitive-local -x ${bow_index} -U ${folder1}/4knext.fasta 2>> ${folder_gral}/RepOthers-ome.log | check_sam ${folder2} | check_mult ${folder2} 100 > ${folder2}/multiple.SAM 
 filter_SAM.sh ${numpro} ${folder2} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 100 &>> ${folder_gral}/RepOthers-ome.log
fi

if [ $startin == "4" ] || [ $startin == "100" ] || [ $startin == "500" ]
then 
 mkdir -p ${folder3}/
 echo " Bowtie2 -k500"
 #Third round
 bowtie2 --seed 22062018  --no-unal --score-min L,0,1.6 -p ${numpro} -k 500 --very-sensitive-local -x ${bow_index} -U ${folder2}/4knext.fastq 2>> ${folder_gral}/RepOthers-ome.log | check_sam ${folder3} | check_mult ${folder3} 500 > ${folder3}/multiple.SAM 
 filter_SAM.sh ${numpro} ${folder3} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 500 &>> ${folder_gral}/RepOthers-ome.log
fi

if [ $startin == "4" ] || [ $startin == "100" ] || [ $startin == "500" ] || [ $startin == "hisat" ]
then
 mkdir -p ${folder4}/
 echo " hisat2"
 #hisat2
 hisat2 --seed 22062018 --un-bz2 ${unali2} -x ${hi_index} -U ${unali} --very-sensitive --novel-splicesite-outfile ${folder_gral}/novel.spli.bed -k 500 -p ${numpro} --no-unal 2>> ${folder_gral}/RepOthers-ome.log | check_sam ${folder4} | check_mult ${folder4} 500 > ${folder4}/multiple.SAM 
 echo -e "\nSpliced Reads:\n" >> ${folder_gral}/summary.txt
 filter_SAM.sh ${numpro} ${folder4} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 500 &>> ${folder_gral}/RepOthers-ome.log
 sort -V -k1,1 -k2,2n ${folder_gral}/novel.spli.bed > ${folder_gral}/splicesites_sorted.bed
 rm ${folder_gral}/novel.spli.bed ${unali}
 echo "Done"
fi


if [ $startin == "4" ] || [ $startin == "100" ] || [ $startin == "500" ] || [ $startin == "hisat" ] || [ $startin == "hisat" ] || [ $startin == "telesc" ]
then
 echo "Mapping Random"
 cat ${folder1}/multiple4random.fastq ${folder2}/multiple4random.fastq ${folder3}/multiple4random.fastq > ${randfastq}
 bowtie2 --seed 22062018 -p ${numpro} -S ${out_rand} --very-sensitive-local --score-min L,0,1.6 --very-sensitive-local -x ${bow_index} -U ${randfastq} &>> ${folder_gral}/RepOthers-ome.log
 samtools view -@ ${numpro} -b ${out_rand} > ${out_rand}.BAM 
 #samtools sort -@ ${numpro} ${out_rand}.BAM -o ${out_rand}_sorted.BAM &>> ${folder_gral}/RepOthers-ome.log
 #samtools merge -f ${out_rand}_all_rand.BAM ${folder1}/DONE_uniq_k4_sorted.BAM ${folder2}/DONE_uniq_k100_sorted.BAM ${folder3}/DONE_uniq_k500_sorted.BAM ${out_rand}_sorted.BAM &>> ${folder_gral}/RepOthers-ome.log
 samtools cat -o ${out_rand}_all_rand.BAM ${folder1}/DONE_uniq_k4.BAM ${folder2}/DONE_uniq_k100.BAM ${folder3}/DONE_uniq_k500.BAM ${out_rand}.BAM &>> ${folder_gral}/RepOthers-ome.log
 samtools sort -@ {numpro} ${out_rand}_all_rand.BAM -o ${out_rand}_all_rand_sorted.BAM &>> ${folder_gral}/RepOthers-ome.log
 rm ${out_rand}_all_rand.BAM ${out_rand}.BAM ${out_rand}
 samtools index ${out_rand}_all_rand_sorted.BAM
 # samtools view -@ ${numpro} ${out_rand}_all_rand.BAM |cut -f 10 > ${folder_gral}/seq4leng.txt 
 # readleng=$(perl -e '{open(SEQ,"$ARGV[0]");$tot=0;$num=0;while($l=<SEQ>){chomp $l;$tot+=length($l);$num++;} $avr=$tot/$num; $avr=int($avr+0.5); print $avr;}' ${folder_gral}/seq4leng.txt)
 # mincov=$(perl -e '{open(SEQ,"$ARGV[0]");$flag=0;while($l=<SEQ>){chomp $l; if($flag==0){$flag=1;$min=length($l);}if(length($l)<$min){$min=length($l);}} print $min;}' ${folder_gral}/seq4leng.txt)
 readleng=$(samtools stats -@ ${numpro} ${out_rand}_all_rand_sorted.BAM | grep -P "^SN\taverage length"| awk '{print $4}')
 mincov=$(samtools stats -@ ${numpro} ${out_rand}_all_rand_sorted.BAM | grep -P "^RL" | sort -V -k2,2n |head -n1|cut -f 2)
 numseq=$(samtools stats -@ ${numpro} ${out_rand}_all_rand_sorted.BAM | grep -P "^SN\tsequences"| cut -f 3)
 FindCoverCutoff.R ${out_rand}_all_rand_sorted.BAM ${readleng} 0.01 &>> ${folder_gral}/RepOthers-ome.log
 cutoff=$(head -n 2 ${out_rand}_all_CoverageCutoff.txt | tail -n 1| cut -f 2) #output de FindCoverCutoff.R
 rm ${randfastq} ${out_rand}*BAM* ${out_rand}_all_CoverageCutoff.txt
 echo "Number of Sequences for RepOthers-ome:" >>${folder_gral}/summary.txt
 echo $numseq >>${folder_gral}/summary.txt
 echo "Average read length:" >>${folder_gral}/summary.txt
 echo $readleng >>${folder_gral}/summary.txt
 echo "Min read length:" >>${folder_gral}/summary.txt
 echo $mincov >>${folder_gral}/summary.txt
 echo "Coverage Cutoff:" >>${folder_gral}/summary.txt
 echo $cutoff >>${folder_gral}/summary.txt
 echo "Done"
 
 echo "Geting possible RepOthers"
 samtools view -H ${folder1}/exons_sorted.BAM|grep -P '(chr\d+|chr\w)\s' | awk '{split($2,a,":"); split($3,b,":"); print a[2],b[2]}' OFS="\t" >${folder_gral}/gen4bedt.txt 
 samtools view -H ${folder1}/exons_sorted.BAM|grep -P '(chr\d+|chr\w)\s' | awk '{split($2,a,":"); split($3,b,":"); print a[2],1,b[2]}' OFS="\t" >${folder_gral}/gen4samt.txt 
 coverage2transcrip.sh ${numpro} ${folder5} ${folder1} ${folder2} ${folder3} ${folder_gral}/gen4samt.txt ${folder_gral}/gen4bedt.txt ${mincov} ${cutoff} &>> ${folder_gral}/RepOthers-ome.log
 echo "Done"

 echo "Telescope"
 network_analysis.sh 10000 75000 ${numpro} ${folder5} ${folder6} ${folder_gral}/gen4samt.txt ${folder_gral}/gen4bedt.txt &>> ${folder_gral}/RepOthers-ome.log 
 rm ${folder_gral}/gen4bedt.txt ${folder_gral}/gen4samt.txt
 echo "Done"
 echo "Final Count"
 awk '{if($5==1)print $1,$2,$3;}' OFS="\t" ${folder5}/regions_filtered_sorted.bed | sort -u > ${folder_gral}/transcripts_solved_telescope.bed
 perl -e '{open(IN,"$ARGV[0]"); while($l=<IN>){ chomp $l; @vec=split("\t",$l); $temp=join("_",@vec);$h{$temp}=1; }close(IN);
  open(FIL,"$ARGV[1]"); while($l=<FIL>){chomp $l; @vec=split("\t",$l); pop @vec; $temp=join("_",@vec); if(!$h{$temp}){print "$l\n";}}}' ${folder_gral}/transcripts_solved_telescope.bed ${folder5}/regions_sorted_coverage_filtered.bed > ${folder_gral}/transcripts_unique.bed
 samtools view -b -L ${folder_gral}/transcripts_unique.bed ${folder5}/ALL.BAM >${folder_gral}/unique.bam 
 bedtools coverage -sorted -mean -a ${folder_gral}/transcripts_solved_telescope.bed -b ${folder6}/result_sorted.bam > ${folder_gral}/transcripts_NOTunique.bed 
 awk -v CUT="$cutoff" '{if (!($4 <= CUT)){ print $0;}}' ${folder_gral}/transcripts_NOTunique.bed > ${folder_gral}/transcripts_NOTunique_filtered.bed
 cat ${folder_gral}/transcripts_NOTunique_filtered.bed ${folder_gral}/transcripts_unique.bed | sort -V -k1,1 -k2,2n > ${folder_gral}/RepOthers.bed
 rm ${folder_gral}/transcripts_NOTunique.bed ${folder_gral}/transcripts_NOTunique_filtered.bed ${folder_gral}/transcripts_unique.bed ${folder_gral}/transcripts_solved_telescope.bed

 samtools merge -f ${folder_gral}/final.bam ${folder6}/result_sorted.bam ${folder_gral}/unique.bam &>> ${folder_gral}/RepOthers-ome.log
 samtools sort -@ ${numpro} ${folder_gral}/final.bam -o ${folder_gral}/RepOthers.bam &>> ${folder_gral}/RepOthers-ome.log
 numseq=$(samtools stats -@ ${numpro} ${folder_gral}/RepOthers.bam | grep -P "^SN\tsequences"| cut -f 3)
 echo "Number of Sequences mapped to RepOthers:" >>${folder_gral}/summary.txt
 echo $numseq >>${folder_gral}/summary.txt
 rm ${folder_gral}/final.bam ${folder_gral}/unique.bam

 awk '{a=$1"_"$2"_"$3; b="gene_id \""a"\"; transcript_id \""a"\"; locus \""a"\";"; print $1,"repeatsome","transcript",$2,$3,".",".",".",b;}' OFS="\t" ${folder_gral}/RepOthers.bed > ${folder_gral}/RepOthers.gtf 
 count_transcripts.R ${folder_gral}/RepOthers.bam ${folder_gral}/RepOthers.gtf transcript gene_id ${folder_gral}/RepOthers_nosplicing FALSE &>> ${folder_gral}/RepOthers-ome.log
 samtools merge -f ${folder_gral}/RepOthers_splicing_temp.bam ${folder4}/DONE_uniq_k500_sorted.BAM ${folder4}/DONE_multiple_k500_sorted.BAM &>> ${folder_gral}/RepOthers-ome.log
 samtools sort -@ ${numpro} ${folder_gral}/RepOthers_splicing_temp.bam -o ${folder_gral}/RepOthers_splicing.bam &>> ${folder_gral}/RepOthers-ome.log
 rm ${folder_gral}/RepOthers_splicing_temp.bam
 count_transcripts.R ${folder_gral}/RepOthers_splicing.bam ${folder_gral}/RepOthers.gtf transcript gene_id ${folder_gral}/RepOthers_splicing TRUE &>> ${folder_gral}/RepOthers-ome.log


 samtools merge -f ${folder_gral}/Exons_temp.bam ${folder1}/exons_sorted.BAM ${folder2}/exons_sorted.BAM ${folder3}/exons_sorted.BAM ${folder4}/exons_sorted.BAM &>> ${folder_gral}/RepOthers-ome.log
 samtools sort -@ ${numpro} ${folder_gral}/Exons_temp.bam -o ${folder_gral}/Exons.bam &>> ${folder_gral}/RepOthers-ome.log
 count_transcripts.R ${folder_gral}/Exons.bam ${exons} exon gene_id ${folder_gral}/Gene FALSE &>> ${folder_gral}/RepOthers-ome.log
 count_transcripts.R ${folder_gral}/Exons.bam ${exons} exon transcript_id ${folder_gral}/Transcript FALSE &>> ${folder_gral}/RepOthers-ome.log

 rm ${folder_gral}/Exons_temp.bam ${folder_gral}/exons_anotation.bed


 rm -r $folder1 $folder2 $folder3 $folder4 $folder5 $folder6

 if [ $paired == "TRUE" ];
 then
  rm ${fastq}
 fi
 echo "Done"
fi

echo "###################   Thanks for using RepOthers-ome   #####################################"



#falta:
#seguros para todos los programas y asi para que pueda ser modular la esto.
#ponerle más y mejor al summary

#Hay que hacer algo con los transcritos? Separar cuales se solucionaron con telescope o no? (info en folder5)
	


#Hay que ponerle seguros por si no tiene instalados diversos paquetes de R o programas tipo bowtie o asi
#bedtools.v.2.25 X
#transcriptR X
#samtools X
#bowtie2 X
#R X
#igraph X
#Rsubread X
#hisat2 X
#parallel

#@article{Tange2011a,
# title = {GNU Parallel - The Command-Line Power Tool},
# author = {O. Tange},
# address = {Frederiksberg, Denmark},
# journal = {;login: The USENIX Magazine},
# month = {Feb},
# number = {1},
# volume = {36},
# url = {http://www.gnu.org/s/parallel},
# year = {2011},
# pages = {42-47}
#}


