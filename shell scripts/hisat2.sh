#!/bin/bash
#SBATCH --job-name=HISAT
#SBATCH --time=0-24:0
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1


for k in 1 2 
  
   do
  
    for z in L005 L006 L007 L008
  
     do
  
      THEMIGBS/softwares/hisat2-2.1.0/hisat2 -p 32 -x GBS/Ref/GRCh38 -1 results/${k}_*_${z}_R1_001.tagged_filter.fastq.gz -2 results/${k}_*_${z}_R2_001.tagged_filter.fastq.gz -S results/SAM/${k}.${z}.sam

	done

  done
