#!/bin/bash

set -e
set -u

FOLDER=$1
NUMPROC=$2
ROUND=$3
GEN4SMT=$4
GEN4BDT=$5
BED=$6


samtools merge -@ ${NUMPROC} -f ${FOLDER}/${ROUND}_temp.bam ${FOLDER}/${ROUND}round_*-updated_sorted.bam
samtools view -@ ${NUMPROC} -b -L ${GEN4SMT} ${FOLDER}/${ROUND}_temp.bam > ${FOLDER}/${ROUND}_2.bam

bedtools intersect -sorted -a ${BED} -b ${FOLDER}/${ROUND}_2.bam -wo -g ${GEN4BDT} > ${FOLDER}/${ROUND}_intersected_reads.bed
sort --parallel ${NUMPROC} -V -u -k 1,3 -k 8,8 ${FOLDER}/${ROUND}_intersected_reads.bed > ${FOLDER}/${ROUND}_intersected_reads_4cov.bed

perl -e '{open(IN,"$ARGV[0]");while(<IN>){@vec=split("\t",$_);$name=$vec[0]."_".$vec[1]."_".$vec[2]; push(@{$h{$vec[7]}}, $name);}close (IN);
 foreach $k (sort keys %h){for ($cont=0;$h{$k}[$cont];$cont++){for($cont2=$cont+1;$h{$k}[$cont2];$cont2++){$out= $h{$k}[$cont]."\t".$h{$k}[$cont2];$done{$out}++;}}} undef %h;
 foreach $k (sort keys %done){print "$k\t$done{$k}\n";} }' ${FOLDER}/${ROUND}_intersected_reads_4cov.bed > ${FOLDER}/${ROUND}_vertex_weight_multiple.txt
cut -f 8 ${FOLDER}/${ROUND}_intersected_reads_4cov.bed | sort| uniq -c | awk '{if($1>1)print $2;}' OFS="\t" >${FOLDER}/${ROUND}_multiple_reads.txt

perl -e '{open(MU, "$ARGV[0]");while($line=<MU>){chomp $line; $mltread{$line}=1;}
 open(IN, "$ARGV[1]");while($line=<IN>){chomp $line; @vec=split("\t",$line); if(!$mltread{$vec[7]}){$name=$vec[0]."_".$vec[1]."_".$vec[2];$h{$name}++;}}close(IN);
 open(VE, "$ARGV[2]"); while($line=<VE>){
  chomp $line; @vec=split("\t",$line); if(!$h{$vec[0]}){$h{$vec[0]}=0;}if(!$h{$vec[1]}){$h{$vec[1]}=0;} $w=$vec[2]+$h{$vec[0]}+$h{$vec[1]};print "$line\t$h{$vec[0]}\t$h{$vec[1]}\t$w\n";}}' ${FOLDER}/${ROUND}_multiple_reads.txt ${FOLDER}/${ROUND}_intersected_reads_4cov.bed ${FOLDER}/${ROUND}_vertex_weight_multiple.txt > ${FOLDER}/vertex_weight.txt

rm ${FOLDER}/${ROUND}round_*-updated_sorted.bam ${FOLDER}/${ROUND}_2.bam ${FOLDER}/${ROUND}_intersected_reads.bed ${FOLDER}/${ROUND}_intersected_reads_4cov.bed ${FOLDER}/${ROUND}_vertex_weight_multiple.txt ${FOLDER}/${ROUND}_multiple_reads.txt


samtools view -@ ${NUMPROC} -F 256 ${FOLDER}/${ROUND}_temp.bam -U ${FOLDER}/${ROUND}_flag256.sam > ${FOLDER}/${ROUND}_flag016.sam
awk '{b=$1","$3","$4; if (!(b in a)){a[b] = $0;} } END { for (i in a) print a[i]}' ${FOLDER}/${ROUND}_flag016.sam > ${FOLDER}/${ROUND}_flag0162.sam
awk '{b=$1","$3","$4; if (!(b in a)){a[b] = $0;} } END { for (i in a) print a[i]}' ${FOLDER}/${ROUND}_flag256.sam > ${FOLDER}/${ROUND}_flag2562.sam
samtools view -H ${FOLDER}/${ROUND}_temp.bam > ${FOLDER}/header.txt


cut -f 1 ${FOLDER}/${ROUND}_flag0162.sam | sort --parallel ${NUMPROC}| uniq -c | awk '{if ($1>1){print $2}}' > ${FOLDER}/${ROUND}_readcorrection.txt
cut -f 1 ${FOLDER}/${ROUND}_flag2562.sam | sort --parallel ${NUMPROC} -u > ${FOLDER}/${ROUND}_readnonuniq.txt
cat ${FOLDER}/${ROUND}_readcorrection.txt >> ${FOLDER}/${ROUND}_readnonuniq.txt


perl -e '{ open(IN,"$ARGV[0]");while($line=<IN>){chomp $line; $h{$line}=1;$hf{$line}=1;}
 open(IN2,"$ARGV[1]"); while(<IN2>){@vec=split("\t",$_);
   if($h{$vec[0]}){
     if($hf{$vec[0]}){
      $hf{$vec[0]}=0;
     }else{
      $vec[1]=$vec[1]+256;
    }} 
   $out=join("\t",@vec); print $out;}}' ${FOLDER}/${ROUND}_readcorrection.txt ${FOLDER}/${ROUND}_flag0162.sam > ${FOLDER}/${ROUND}_modified2.sam


perl -e '{ open(IN,"$ARGV[0]");while($line=<IN>){chomp $line; $h{$line}=1;}
 open(IN2,"$ARGV[1]"); while(<IN2>){@vec=split("\t",$_);
   if(!$h{$vec[0]}){print $_;}
   }}' ${FOLDER}/${ROUND}_readnonuniq.txt ${FOLDER}/${ROUND}_flag0162.sam> ${FOLDER}/${ROUND}_uniqreads.sam


ROUND2=$(($ROUND+1))

cat ${FOLDER}/${ROUND}_modified2.sam ${FOLDER}/${ROUND}_flag2562.sam | awk '{b=$1","$3","$4; if (!(b in a)){a[b] = $0;} } END { for (i in a) print a[i]}'| sort -V -k1,1 -k2,2n  > ${FOLDER}/${ROUND2}.SAM


parallel --pipepart -k --block -1 -j ${NUMPROC} -a ${FOLDER}/${ROUND2}.SAM -q awk '{if ($13~/XS/){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21;}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20;}}' OFS="\t" > ${FOLDER}/${ROUND}_modified2.sam
cat ${FOLDER}/header.txt ${FOLDER}/${ROUND}_modified2.sam > ${FOLDER}/${ROUND2}.SAM
rm  ${FOLDER}/${ROUND}_readnonuniq.txt ${FOLDER}/${ROUND}_modified2.sam ${FOLDER}/${ROUND}_readcorrection.txt ${FOLDER}/${ROUND}_flag016.sam ${FOLDER}/${ROUND}_flag256.sam ${FOLDER}/${ROUND}_flag2562.sam ${FOLDER}/${ROUND}_flag0162.sam


