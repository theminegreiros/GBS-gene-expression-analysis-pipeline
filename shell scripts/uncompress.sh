#!/bin/bash
for k in 13
 do
  for n in {5..8}
   do
   for z in 1 2 
    do
    # gzip -d /data3/theminegreiros/GBS/NPAD/decompressed/${k}/${k}_*_L00${n}_R${z}_001.fastq.gz
      pigz -p 16 -d /data3/theminegreiros/GBS/NPAD/decompressed/${k}_*_L00${n}_R${z}_001.fastq.gz
   done
  done
 done
