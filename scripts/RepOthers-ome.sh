#!/bin/bash

set -e
set -u
numpro=1
folder_gral=./
bow_index=""
rtRNA_MChr=""
exons=""
hi_index=""
fastq1=""
fastq2=""
fastq=""

while getopts 'p:b:r:e:o:i:1:2:U:h' OPTION;do
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
 ?)
  echo "script usage: $(basename $0) [-l] [-h] [-a somevalue]" >&2 #mensaje de ayuda para correr comando en caso de que haya un parametro que no venga al caso
  exit 1
  ;;

 esac
done

if [ -z "$bow_index" ] || [ -z "$rtRNA_MChr" ] || [ -z "$exons" ] || [ -z "$hi_index" ] 
then
 echo "Mensaje de error 2" && exit 1
fi

if [ -z "$fastq1" ] || [ -z "$fastq2" ]
then
 if [ -z "$fastq" ]
 then
  echo "Mensaje de error 3" && exit 1
 else
  paired=FALSE	
 fi
else
 paired=TRUE
fi



mkdir -p ${folder_gral}/
if [ $paired == "TRUE" ];
then
 fastq=${folder_gral}/test_joined.fastq
 >${fastq}
 cat ${fastq1} | parallel --tmpdir ${folder_gral} --block 500M -j ${numpro} --pipe -L 4000000 --cat perl\ -e\ \'\{open\(IN,\$ARGV\[0\]\)\;\ while\(\<IN\>\)\{@vec\=split\(\"\ \",\$_\)\;\$l\=\$vec\[0\].\"_\".\$ARGV\[1\].\"\\n\"\;\ \$a\=\<IN\>\;\$b\=\<IN\>\;\$c\=\<IN\>\;print\ \"\$l\$a+\\n\$c\"\;\}\}\'\ \{\}\ 1 >>${fastq}
 cat ${fastq2} | parallel --tmpdir ${folder_gral} --block 500M -j ${numpro} --pipe -L 4000000 --cat perl\ -e\ \'\{open\(IN,\$ARGV\[0\]\)\;\ while\(\<IN\>\)\{@vec\=split\(\"\ \",\$_\)\;\$l\=\$vec\[0\].\"_\".\$ARGV\[1\].\"\\n\"\;\ \$a\=\<IN\>\;\$b\=\<IN\>\;\$c\=\<IN\>\;print\ \"\$l\$a+\\n\$c\"\;\}\}\'\ \{\}\ 2 >>${fastq}
fi





unali=${folder_gral}/test_unalign.bz2
out_1=${folder_gral}/k4.SAM
out_2=${folder_gral}/k100.SAM
out_3=${folder_gral}/k500.SAM
out_4=${folder_gral}/hisat_k500.SAM
folder1=${folder_gral}/mapping_4
folder2=${folder_gral}/mapping_100
folder3=${folder_gral}/mapping_500
folder4=${folder_gral}/mapping_splicesites
folder5=${folder_gral}/transcripts
folder6=${folder_gral}/network
randfastq=${folder_gral}/for_random.fq
out_rand=${folder_gral}/random




parallel --pipepart -a ${exons} -j ${numpro} --block -1 -q awk '{if($3 == "exon")print $1,$4,$5;}' OFS="\t" ${exons} | sort --parallel ${numpro} -V -k1,1 -k2,2n -u >${folder_gral}/exons_anotation.bed
>${folder_gral}/summary.txt
>${folder_gral}/RepOthers-ome.log

echo "Mapping"
echo " Bowtie2 -k4"
#First Bowtie2 with a small number of possible alignments
bowtie2 --un-bz2 ${unali} --no-unal --score-min L,0,1.6 -p ${numpro} -k 4 -S ${out_1} --very-sensitive-local -x ${bow_index} -U ${fastq} &>> ${folder_gral}/RepOthers-ome.log
filter_SAM.sh ${numpro} ${folder1} ${out_1} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 4 1>> ${folder_gral}/RepOthers-ome.log
#rm ${out_1}
echo " Bowtie2 -k100"
#Second round
bowtie2 --no-unal --score-min L,0,1.6 -p ${numpro} -k 100 -S ${out_2} --very-sensitive-local -x ${bow_index} -U ${folder1}/4knext.fastq &>> ${folder_gral}/RepOthers-ome.log
filter_SAM.sh ${numpro} ${folder2} ${out_2} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 100 1>> ${folder_gral}/RepOthers-ome.log
#rm ${out_2}
echo " Bowtie2 -k500"
#Third round
bowtie2 --no-unal --score-min L,0,1.6 -p ${numpro} -k 500 -S ${out_3} --very-sensitive-local -x ${bow_index} -U ${folder2}/4knext.fastq &>> ${folder_gral}/RepOthers-ome.log
filter_SAM.sh ${numpro} ${folder3} ${out_3} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 500 1>> ${folder_gral}/RepOthers-ome.log
#rm ${out_3}
echo " hisat2"
#hisat2
hisat2 -x ${hi_index} -U ${unali} --very-sensitive --novel-splicesite-outfile ${folder_gral}/novel.spli.bed -k 500 -p ${numpro} --no-unal -S ${out_4} &>> ${folder_gral}/RepOthers-ome.log
echo -e "Spliced Reads:\n\n" >> ${folder_gral}/summary.txt
filter_SAM.sh ${numpro} ${folder4} ${out_4} ${folder_gral}/summary.txt ${folder_gral}/exons_anotation.bed ${rtRNA_MChr} 500 1>> ${folder_gral}/RepOthers-ome.log
sort -V -k1,1 -k2,2n ${folder_gral}/novel.spli.bed > ${folder_gral}/splicesites_sorted.bed
#rm ${out_4}
rm ${folder_gral}/novel.spli.bed ${unali}
echo "Done"


echo "Mapping Random"
cat ${folder1}/multiple4random.fastq ${folder2}/multiple4random.fastq ${folder3}/multiple4random.fastq > ${randfastq}
bowtie2 -p ${numpro} -S ${out_rand} --very-sensitive-local --score-min L,0,1.6 --very-sensitive-local -x ${bow_index} -U ${randfastq} &>> ${folder_gral}/RepOthers-ome.log
samtools view -@ ${numpro} -b ${out_rand} > ${out_rand}.BAM 
samtools sort -@ ${numpro} ${out_rand}.BAM -o ${out_rand}_sorted.BAM &>> ${folder_gral}/RepOthers-ome.log
rm ${out_rand}.BAM ${out_rand}
samtools merge -f ${out_rand}_all_rand.BAM ${folder1}/DONE_uniq_k4_sorted.BAM ${folder2}/DONE_uniq_k100_sorted.BAM ${folder3}/DONE_uniq_k500_sorted.BAM ${out_rand}_sorted.BAM &>> ${folder_gral}/RepOthers-ome.log
samtools index ${out_rand}_all_rand.BAM
samtools view -@ ${numpro} ${out_rand}_all_rand.BAM |cut -f 10 > ${folder_gral}/seq4leng.txt 
readleng=$(perl -e '{open(SEQ,"$ARGV[0]");$tot=0;$num=0;while($l=<SEQ>){chomp $l;$tot+=length($l);$num++;} $avr=$tot/$num; $avr=int($avr+0.5); print $avr;}' ${folder_gral}/seq4leng.txt)
mincov=$(perl -e '{open(SEQ,"$ARGV[0]");$flag=0;while($l=<SEQ>){chomp $l; if($flag==0){$flag=1;$min=length($l);}if(length($l)<$min){$min=length($l);}} print $min;}' ${folder_gral}/seq4leng.txt)
FindCoverCutoff.R ${out_rand}_all_rand.BAM ${readleng} 0.01 &>> ${folder_gral}/RepOthers-ome.log
cutoff=$(head -n 2 ${out_rand}_CoverageCutoff.txt | tail -n 1| cut -f 2) #output de FindCoverCutoff.R
rm ${folder_gral}/seq4leng.txt ${randfastq} ${out_rand}*BAM* ${out_rand}_CoverageCutoff.txt
echo "Done"


samtools view -H ${folder1}/exons_sorted.BAM|grep -P '(chr\d+|chr\w)\s' | awk '{split($2,a,":"); split($3,b,":"); print a[2],b[2]}' OFS="\t" >${folder_gral}/gen4bedt.txt 
samtools view -H ${folder1}/exons_sorted.BAM|grep -P '(chr\d+|chr\w)\s' | awk '{split($2,a,":"); split($3,b,":"); print a[2],1,b[2]}' OFS="\t" >${folder_gral}/gen4samt.txt 

echo "Geting possible RepOthers"
coverage2transcrip.sh ${numpro} ${folder5} ${folder1} ${folder2} ${folder3} ${folder_gral}/gen4samt.txt ${folder_gral}/gen4bedt.txt ${mincov} ${cutoff} &>> ${folder_gral}/RepOthers-ome.log
echo "Done"

echo "Telescope"
network_analysis.sh 1500 ${numpro} ${folder5} ${folder6} ${folder_gral}/gen4samt.txt ${folder_gral}/gen4bedt.txt &>> ${folder_gral}/RepOthers-ome.log 
rm ${folder_gral}/gen4bedt.txt ${folder_gral}/gen4samt.txt
echo "Done"

echo "Final Count"
awk '{if($5==1)print $1,$2,$3;}' OFS="\t" ${folder5}/regions_filtered_sorted.bed | sort -u > ${folder_gral}/transcripts_solved_telescope.bed
perl -e '{open(IN,"$ARGV[0]"); while($l=<IN>){ chomp $l; @vec=split("\t",$l); $temp=join("_",@vec);$h{$temp}=1; }close(IN);
 open(FIL,"$ARGV[1]"); while($l=<FIL>){chomp $l; @vec=split("\t",$l); pop @vec; $temp=join("_",@vec); if(!$h{$temp}){print "$l\n";}}}' ${folder_gral}/transcripts_solved_telescope.bed ${folder5}/regions_sorted_coverage_filtered.bed > ${folder_gral}/transcripts_unique.bed
samtools view -b -L ${folder_gral}/transcripts_unique.bed ${folder5}/ALL.BAM >${folder_gral}/unique.bam 
bedtools coverage -mean -a ${folder_gral}/transcripts_solved_telescope.bed -b ${folder6}/result_sorted.bam > ${folder_gral}/transcripts_NOTunique.bed 
awk -v CUT="$cutoff" '{if (!($4 <= CUT)){ print $0;}}' ${folder_gral}/transcripts_NOTunique.bed > ${folder_gral}/transcripts_NOTunique_filtered.bed
cat ${folder_gral}/transcripts_NOTunique_filtered.bed ${folder_gral}/transcripts_unique.bed | sort -V -k1,1 -k2,2n > ${folder_gral}/RepOthers.bed
rm ${folder_gral}/transcripts_NOTunique.bed ${folder_gral}/transcripts_NOTunique_filtered.bed ${folder_gral}/transcripts_unique.bed ${folder_gral}/transcripts_solved_telescope.bed

samtools merge -f ${folder_gral}/final.bam ${folder6}/result_sorted.bam ${folder_gral}/unique.bam &>> ${folder_gral}/RepOthers-ome.log
samtools sort -@ ${numpro} ${folder_gral}/final.bam -o ${folder_gral}/RepOthers.bam &>> ${folder_gral}/RepOthers-ome.log
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
echo "###################   Thanks for using RepOthers-ome   #####################################"



#falta:
#hacer correr esto
#ponerle parámetros
#seguros para todos los programas y asi para que pueda ser modular la esto.
#ponerle más y mejor al summary
#AMI A LO QUE TENEMOS!!!!

#Hay que hacer algo con los transcritos? Separar cuales se solucionaron con telescope o no? (info en folder5)



#Hay que ponerle seguros por si no tiene instalados diversos paquetes de R o programas tipo bowtie o asi
#bedtools.v.2.25
#transcriptR
#samtools
#bowtie2
#R
#igraph
#Rsubread
#hisat2
#parallel


#O. Tange (2011): GNU Parallel - The Command-Line Power Tool,
#  ;login: The USENIX Magazine, February 2011:42-47.

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


