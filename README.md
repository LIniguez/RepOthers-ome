# RepOthers-ome
A methodology that finds RepOthers. It is based on bowtie2 and telescope

The sort of bowtie2 index is very very important since from it depends the corrct use of RepOthers-ome.
It should be on the way chr1,chr2,chr3... (no underscores for primary chromosomes) and they shoud be sorted as version mode (sort -V).
The hisat2 index needs to have the same chromosomal names but not the same order. 

It is very important the annotation file for the exons, chromosomes should be equal in the indexes. It is recomended only chromosomal
annotations. 