#!/bin/bash
#SBATCH --job-name=samtobam
#SBATCH --time=10-0:0
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

for k in 1

  do
    for z in L005 L006 L007 L008
     do
      softwares/samtools-1.3.1/samtools view -u SAM/${k}.${z}.sam | softwares/samtools-1.3.1/samtools sort -@ 32 > BAM/${k}.${z}.bam
                                
     done
  done
