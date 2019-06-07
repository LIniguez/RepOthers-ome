#!/bin/bash

set -e
set -u

FOLDOUT=$1
NPRO=$2
I=$3
GEN4ST=$4
GEN4BT=$5
BED=$6
THRE=$7

samtools merge -@ ${NPRO} -f ${FOLDOUT}/${I}_temp.bam ${FOLDOUT}/${I}round_*-updated_sorted.bam
samtools view -H ${FOLDOUT}/${I}_temp.bam | grep -P '@SQ'| awk '{split($2,a,":"); split($3,b,":"); print a[2],b[2]}' OFS="\t" >${FOLDOUT}/header_mod.txt
#samtools view -H ${FOLDOUT}/${I}_temp.bam >${FOLDOUT}/header_mod.txt
samtools view -@ ${NPRO} -L ${GEN4ST} ${FOLDOUT}/${I}_temp.bam | awk THRE=${THRE}'{ for (i=1; i<=NF; ++i) { if ($i ~ "XP:i") {split($i,a,":");if(a[3]>THRE){ if($0 !~ "ZF:Z:__no_feature"){print $0;next;}}}}}' OFS="\t" > ${FOLDOUT}/${I}_telescope_res.sam #THRE is the threshold for telescope probability of assignment. 
samtools view -@ ${NPRO} -t ${FOLDOUT}/header_mod.txt -b ${FOLDOUT}/${I}_telescope_res.sam > ${FOLDOUT}/${I}_2.bam  #remove sequence not present in gen4smt regions

bedtools intersect -sorted -a ${BED} -b ${FOLDOUT}/${I}_2.bam -wo -g ${GEN4BT} | sort --parallel ${NPRO} -V -u -k 1,3 -k 8,8 > ${FOLDOUT}/${I}_intersected_reads_4cov.bed

perl -e '{open(IN,"$ARGV[0]");while(<IN>){@vec=split("\t",$_);$name=$vec[0]."_".$vec[1]."_".$vec[2]; push(@{$h{$vec[7]}}, $name);}close (IN);
 foreach $k (sort keys %h){for ($cont=0;$h{$k}[$cont];$cont++){for($cont2=$cont+1;$h{$k}[$cont2];$cont2++){$out= $h{$k}[$cont]."\t".$h{$k}[$cont2];$done{$out}++;}}} undef %h;
 foreach $k (sort keys %done){print "$k\t$done{$k}\n";} }' ${FOLDOUT}/${I}_intersected_reads_4cov.bed > ${FOLDOUT}/${I}_vertex_weight_multiple.txt
cut -f 8 ${FOLDOUT}/${I}_intersected_reads_4cov.bed | sort| uniq -c | awk '{if($1>1)print $2;}' OFS="\t" >${FOLDOUT}/${I}_multiple_reads.txt

perl -e '{open(MU, "$ARGV[0]");while($line=<MU>){chomp $line; $mltread{$line}=1;}
 open(IN, "$ARGV[1]");while($line=<IN>){chomp $line; @vec=split("\t",$line); if(!$mltread{$vec[7]}){$name=$vec[0]."_".$vec[1]."_".$vec[2];$h{$name}++;}}close(IN);
 open(VE, "$ARGV[2]"); while($line=<VE>){
  chomp $line; @vec=split("\t",$line); if(!$h{$vec[0]}){$h{$vec[0]}=0;}if(!$h{$vec[1]}){$h{$vec[1]}=0;} $w=$vec[2]+$h{$vec[0]}+$h{$vec[1]};print "$line\t$h{$vec[0]}\t$h{$vec[1]}\t$w\n";}}' ${FOLDOUT}/${I}_multiple_reads.txt ${FOLDOUT}/${I}_intersected_reads_4cov.bed ${FOLDOUT}/${I}_vertex_weight_multiple.txt > ${FOLDOUT}/vertex_weight.txt

rm ${FOLDOUT}/${I}round_*-updated_sorted.bam   ${FOLDOUT}/${I}_intersected_reads_4cov.bed ${FOLDOUT}/${I}_vertex_weight_multiple.txt ${FOLDOUT}/${I}_multiple_reads.txt



samtools view -@ ${NPRO} -F 256 ${FOLDOUT}/${I}_2.bam -U ${FOLDOUT}/${I}_flag256.sam > ${FOLDOUT}/${I}_flag016.sam
awk '{b=$1","$3","$4; if (!(b in a)){a[b] = $0;} } END { for (i in a) print a[i]}' ${FOLDOUT}/${I}_flag016.sam > ${FOLDOUT}/${I}_flag0162.sam &
awk '{b=$1","$3","$4; if (!(b in a)){a[b] = $0;} } END { for (i in a) print a[i]}' ${FOLDOUT}/${I}_flag256.sam > ${FOLDOUT}/${I}_flag2562.sam &
wait


cut -f 1 ${FOLDOUT}/${I}_flag0162.sam | sort --parallel ${NPRO}| uniq -c | awk '{if ($1>1){print $2}}' > ${FOLDOUT}/${I}_readcorrection.txt &
cut -f 1 ${FOLDOUT}/${I}_flag2562.sam | sort --parallel ${NPRO} -u > ${FOLDOUT}/${I}_readnonuniq.txt &
wait
cat ${FOLDOUT}/${I}_readcorrection.txt >> ${FOLDOUT}/${I}_readnonuniq.txt


perl -e '{ open(IN,"$ARGV[0]");while($line=<IN>){chomp $line; $h{$line}=1;$hf{$line}=1;}
 open(IN2,"$ARGV[1]"); while(<IN2>){@vec=split("\t",$_);
   if($h{$vec[0]}){
     if($hf{$vec[0]}){
      $hf{$vec[0]}=0;
     }else{
      $vec[1]=$vec[1]+256;
    }}
   $out=join("\t",@vec); print $out;}}' ${FOLDOUT}/${I}_readcorrection.txt ${FOLDOUT}/${I}_flag0162.sam > ${FOLDOUT}/${I}_modified2.sam


perl -e '{ open(IN,"$ARGV[0]");while($line=<IN>){chomp $line; $h{$line}=1;}
 open(IN2,"$ARGV[1]"); while(<IN2>){@vec=split("\t",$_);
   if(!$h{$vec[0]}){print $_;}
   }}' ${FOLDOUT}/${I}_readnonuniq.txt ${FOLDOUT}/${I}_flag0162.sam> ${FOLDOUT}/${I}_uniqreads.sam


I2=$(($I+1))

cat ${FOLDOUT}/${I}_modified2.sam ${FOLDOUT}/${I}_flag2562.sam | awk '{b=$1","$3","$4; if (!(b in a)){a[b] = $0;} } END { for (i in a) print a[i]}'| sort -T ${FOLDOUT} --parallel ${NPRO} -V -k1,1 -k2,2n  > ${FOLDOUT}/${I2}.SAM


awk '{if ($13~/XS/){print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21;}else{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20;}}' OFS="\t" ${FOLDOUT}/${I2}.SAM > ${FOLDOUT}/${I}_modified2.sam
samtools view -@ ${NPRO} -b ${FOLDOUT}/${I}_modified2.sam -t ${FOLDOUT}/header_mod.txt > ${FOLDOUT}/${I2}_temp.BAM
rm  ${FOLDOUT}/${I}_temp.bam ${FOLDOUT}/${I}_readnonuniq.txt ${FOLDOUT}/header_mod.txt ${FOLDOUT}/${I}_2.bam ${FOLDOUT}/${I}_modified2.sam ${FOLDOUT}/${I}_readcorrection.txt ${FOLDOUT}/${I}_flag016.sam ${FOLDOUT}/${I}_flag256.sam ${FOLDOUT}/${I}_flag2562.sam ${FOLDOUT}/${I}_flag0162.sam ${FOLDOUT}/${I2}.SAM
