#!/bin/bash
for k in {15..24}
 do
  for n in {5..8}
   do
    /data/home/theminegreiros/sortmerna-3.0.3/scripts/unmerge-paired-reads.sh ${k}_L00${n}_reads_non_rRNA.fastq  ${k}_L00${n}_reads_non_rRNA_1.fastq  ${k}_L00${n}_reads_non_rRNA_2.fastq
    gzip ${k}_L00${n}_reads_non_rRNA.fastq 
    gzip ${k}_L00${n}_reads_non_rRNA_1.fastq
    gzip ${k}_L00${n}_reads_non_rRNA_2.fastq
    
   done
  done
