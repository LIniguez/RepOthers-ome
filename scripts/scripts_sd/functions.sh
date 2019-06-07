#!/bin/sh

check_sam(){
 >${1}unique.SAM; >${1}best_alscor.txt; >${1}header.txt
 awk -v FOLD="$1" '{
 if($1~/^@/) {print $0 >> FOLD"header.txt";}                              #Remove header sequencs
 else if(($2-256)<0){                                                     #check if the alignment is the principal
  if($13~/XS/){                                                           #if the alignment has XS it means it has multiple positions in Bowtie2
    split($12,a,":");print $0 ; print $1,a[3] >> FOLD"best_alscor.txt";   #print it to multiple.SAM and retain the info of the alignement score
  }else if($13~/ZS/){                                                     #if the alignment has ZS it means it has multiple positions in Bowtie2
     split($12,a,":"); print $0 ; print $1,a[3] >> FOLD"best_alscor.txt"; #print it to multiple.SAM and retain the info of the alignement score
  }else{print $0 >> FOLD"unique.SAM";}}
  else{ print $0 ;}}' OFS="\t"
}
check_mult(){
 FOLD=$1; MAX=$2
 perl -pe '{@vec=split("\t",$_);$h{$vec[0]}++;END{open(OUT, ">$OUT");
 foreach $read(keys %h){if( $h{$read} < $MAX ){ print OUT "$read\n";}}}}' -s -- -OUT=${FOLD}multreads_done.txt -MAX=${MAX}
}
logo(){
  printf '\n'
  printf '#                                                   *\n'
  printf '##                  *         888888888888888888888888888     *\n'
  printf '###       *                   88      88    88 88   88  88\n'
  printf '####                       *   8888   88   88   88  88888\n'
  printf '############       *              88  88  888888888 88  88\n'
  printf '##########                  8888888   88  88     88 88   88\n'
  printf '########     *                       *\n'
  printf '######                 *    88888888888  88888888888888888\n'
  printf '#####            *          88   88  88  88  88  *  88\n'
  printf '######                      88    88 88  88   8888  88\n'
  printf '  ######     *      *       88   88  88  88      88 88\n'
  printf '     ####                   88888     8888   88888  88     *\n'
  printf '        ##                       *\n'
  printf '\n'
}
