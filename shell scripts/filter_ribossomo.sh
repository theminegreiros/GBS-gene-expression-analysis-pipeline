#!/bin/bash
for k in {7..12}
 do
  for n in {5..8}
   do
    sortmerna --ref /data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/silva-bac-16s-id90.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/silva-bac-16s-db:/data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/silva-bac-23s-id98.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/silva-bac-23s-db:/data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/silva-arc-16s-id95.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/silva-arc-16s-db:/data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/silva-arc-23s-id98.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/silva-arc-23s-db:/data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/silva-euk-18s-id95.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/silva-euk-18s-db:/data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/silva-euk-28s-id98.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/silva-euk-28s:/data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/rfam-5s-database-id98.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/rfam-5s-db:/data/home/theminegreiros/sortmerna-3.0.3/rRNA_databases/rfam-5.8s-database-id98.fasta,/data/home/theminegreiros/sortmerna-3.0.3/index/rfam-5.8s-db --reads /data3/theminegreiros/GBS/NPAD/decompressed/${k}/${k}_L00${n}_*.fastq --num_alignments 1 --fastx --aligned ${k}_L00${n}_reads_rRNA --other ${k}_L00${n}_reads_non_rRNA --log -a 16 -m 8192 --paired_in -v
    done
  done

