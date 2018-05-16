#!/bin/bash
FOLDER=$1
ROUND=$2


for i in $(find ${FOLDER}*_${ROUND}.txt)  ##### change
do
SRR=$(echo $i| rev | cut -f2- -d "."| rev)
sed -i 's/\r$//' ${i}
awk '{split($2,a,"_");if(a[2]>0){b="gene_id \""$2"\"; transcript_id \""$2"\"; locus \""$2"\";"; print a[1],"repeatsome","transcript",a[2],a[3],".",".",".",b;}}' OFS="\t" ${SRR}.txt > ${SRR}.gtf
rm ${i}
done



