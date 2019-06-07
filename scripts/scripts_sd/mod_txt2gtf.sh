#!/bin/bash

set -e
set -u


FOLDOUT=$1
I=$2


for i in $(find ${FOLDOUT}*_${I}.txt)  ##### change
do
NAME=$(echo $i| rev | cut -f2- -d "."| rev)
sed -i 's/\r$//' ${i}
awk '{split($2,a,"_");if(a[2]>0){b="gene_id \""$2"\"; transcript_id \""$2"\"; locus \""$2"\";"; print a[1],"SDtemp","transcript",a[2],a[3],".",".",".",b;}}' OFS="\t" ${NAME}.txt > ${NAME}.gtf
rm ${i}
done
