#!/bin/sh

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

while getopts 'pb:r:e:oi:1:2:U:h' OPTION;do
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

# numpro=$(nproc)   -p
# paired=TRUE -P #este queda implÃ­cito si hay -1 y -2 o solo -U
# bow_index=/home/ubuntu/references/hg38/bowtie2/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index -b
# rtRNA_MChr=/home/ubuntu/tests/hg38_rtM_RNA.bed -r
# exons=/home/ubuntu/tests/gencode.v28.annotation.gtf -e
# folder_gral=/home/ubuntu/tests -o
# hi_index=/home/ubuntu/efs/references/hisat2/hg38/genome -i
# help -h
# _1.fastq -1
# _2.fastq -2
# fastq -U

