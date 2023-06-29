#!/bin/bash

for i in `ls`; 
do 
zcat $i/*1.fq.gz >> ../combined/$i.1.fq; 
echo $i.1 ;
zcat $i/*2.fq.gz >> ../combined/$i.2.fq;
echo $i.2 ; 
STAR --genomeDir ../NCBI_hg38_star_149bp_overhang --runThreadN 24 --readFilesIn ../combined/$i.1.fq ../combined/$i.2.fq --outFileNamePrefix ../result/$i/$i_ --sjdbOverhang 149 --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts;
echo $i
rm ../combined/$i.1.fq ../combined/$i.2.fq;
rm $i/*1.fq.gz $i/*2.fq.gz;
done
