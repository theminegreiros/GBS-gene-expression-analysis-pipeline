#!/bin/bash
for k in {1..12}
 do
  for n in {5..8}
   do
    #for z in 1 2
     #do 
      # gzip ${k}_*_L00${n}_R${z}_001.fastq
        gzip ${k}/${k}_L00${n}_interleaved.* 
     #done
    done
 done

