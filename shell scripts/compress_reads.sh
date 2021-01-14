#!/bin/bash
for k in {1..14}
 do
  for n in {5..8}
   do
     gzip ${k}_L00${n}_reads_non_rRNA.fastq
     gzip ${k}_L00${n}_reads_non_rRNA_1.fastq
     gzip ${k}_L00${n}_reads_non_rRNA_2.fastq
     gzip ${k}_L00${n}_reads_rRNA.fastq	  	

   done
 done
