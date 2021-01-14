# GBS gene expression analysis pipeline (RNA-Seq)
Gene expression analysis of peripheral blood cells (PBMC) of Guillian-Barré Syndrome Patients.
## Pipeline steps
## Linux environment
### Softwares and tools descriptions
* [Fastp](https://github.com/OpenGene/fastp) v.0.20.0:  FASTQ data pre-processing tool. The algorithm has functions for quality control, trimming of adapters, filtering by quality, and read pruning.
* [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) v.0.11.9: A quality control tool for high throughput sequence data.
* [FeatureCounts](http://bioinf.wehi.edu.au/featureCounts/) v.1.6.5: A highly efficient general-purpose read summarization program that counts mapped reads for genomic features such as genes, exons, promoter, gene bodies, genomic bins and chromosomal locations.
* [Hisat2](http://daehwankimlab.github.io/hisat2/) v.2.2.0: A fast and sensitive alignment program for mapping next-generation sequencing reads (both DNA and RNA) to a population of human genomes as well as to a single reference genome.
* [SAMtools](http://www.htslib.org/doc/samtools.html) v.1.10: A set of utilities that manipulate alignments in the BAM format. It imports from and exports to the SAM (Sequence Alignment/Map) format, does sorting, merging and indexing, and allows to retrieve reads in any regions swiftly.
* [SortMeRNA](https://bioinfo.lifl.fr/RNA/sortmerna/) v.2.1: A program tool for filtering, mapping and OTU-picking NGS reads in metatranscriptomic and metagenomic data.
### 1. Check quality of transcriptome libraries
  * FastQC v.0.11.9
### 2. Filter Ribossomal RNA
  * SortMeRNA v.2.1
### 3. Trimming, adapters and low quality reads removal  
  * Fastp v.0.20.0
### 4. Check quality of filtered transcriptome libraries  
  * FastQC v.0.11.9
### 5. Mapping reads to the Human genome (GRCh38)
  * Hisat2 v.2.2.0
### 6. Converting Sequence Alignment/Map (SAM) to Binary Alignment Map (BAM)
  * SAMtools v.1.10 
### 7. Count mapped reads
  * FeatureCounts v.0.11.9
## R environment
### R packages descriptions
* [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html): Provide access to a diverse set of data and enables a wide range of powerful online queries from gene annotation to database mining.
* [clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) : This package implements methods to analyze and visualize functional profiles (GO and KEGG) of gene and gene clusters.
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html): Estimate variance-mean dependence in count data from high-throughput sequencing assays and test for differential expression based on a model using the negative binomial distribution.
* [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html): Genome wide annotation for Human, primarily based on mapping using Entrez Gene identifiers.
* [pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12/topics/pheatmap): A function to draw clustered heatmaps where one has better control over some graphical parameters such as cell size, etc.
* [genefilter](https://bioconductor.org/packages/release/bioc/html/genefilter.html): Methods for filtering genes from high-throughput experiments.
* [ggplot2](https://github.com/tidyverse/ggplot2): A system for declaratively creating graphics, based on The Grammar of Graphics.
* [RColorBrewer](https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2/topics/RColorBrewer): Creates nice looking color palettes especially for thematic maps.

### 8. Load, check and preprocess the count mapped dataframe
### 9. Differential Gene Expression
* DESeq2
