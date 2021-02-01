#!/bin/bash
#SBATCH --job-name=COUNTS
#SBATCH --time=10-0:0
#SBATCH --mail-user=theminegreiros@yahoo.com.br
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=32
#SBATCH --nodes=1

THEMIGBS/softwares/subread-1.6.2-Linux-x86_64/bin/featureCounts -T 32 -p -t gene -g gene_id -a GBS/Ref/Homo_sapiens.GRCh38.92.gff3 -o counts/counts_npad.txt results/BAM/*.bam
