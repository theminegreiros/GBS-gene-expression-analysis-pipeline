# GBS gene expression analysis pipeline
Gene expression analysis of peripheral blood cells (PBMC) of Guillian-Barré Syndrome Patients.
## Pipeline steps
### 1. Check quality of transcriptome libraries
  * FastQC v.0.11.9
### 2. Filter Ribossomal RNA
  * SortMeRNA v.2.1
### 3. Trimming, adapters and low quality reads removal  
  * Fastp v.0.20.0
### 4. Mapping reads to the human genome
  * Hisat2 v.2.2.0
