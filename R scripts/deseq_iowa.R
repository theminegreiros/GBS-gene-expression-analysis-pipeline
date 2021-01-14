# setwd("/home/themi/IOWA_notebook/counts")
setwd("D:/IOWA_notebook/counts")
library(DESeq2)
library("ggplot2")
#library("vsn")
library("pheatmap")
library("RColorBrewer")
library("org.Hs.eg.db")
library("clusterProfiler")
counts_hungria <- read.table("counts_non_rRNA_iowa.txt", header = T, stringsAsFactors = F, sep = "\t")
countshungria = counts_hungria[,c(1,7:ncol(counts_hungria))]
counts_npad <- read.table("counts_non_rRNA_NPAD.txt", header = T, stringsAsFactors = F, sep = "\t")
counts = merge(counts_npad,countshungria, by = "Geneid")
cts <- counts[, c(1, 7:ncol(counts))]
cts <- cts[,-c(2:9, 18:21, 26:29,38:41, 42:45, 78:93 )]

cond = c(#rep("REC60", 4),#10##Demyelinating
  #rep("GBS", 4),#11##Demyelinating
  # rep("T2", 4),#12##Miller_Fisher##120
  # rep("T1", 4),#13##Axonal
  # #rep("REC60", 4),#14##Axonal
  # rep("T1", 4),#15##Demyelinating
  # #rep("GBS", 4),#16##Miller_Fisher
  # rep("T2", 4),#17##Miller_Fisher##60
  # rep("T1", 4),#18##Axonal
  # #rep("REC120", 4),#19##Axonal
  # #rep("GBS", 4),#1##Demyelinating
  # rep("T2", 4),#20##Demyelinating##180
  # rep("T1", 4),#21##Demyelinating
  # rep("T2", 4),#22##Demyelinating##240
  # rep("T2", 4),#23##Demyelinating##60
  # rep("T1", 4),#24##Demyelinating
  # rep("T2", 4),#2##Demyelinating##60
  # rep("T1", 4),#3##Demyelinating
  # rep("T2", 4),#4##Demyelinating##120
  # #rep("GBS", 4),#5##Axonal
  # #rep("REC60", 4),#6##Axonal
  # #rep("GBS", 4),#7##Axonal
  # #rep("REC90", 4),#8##Axonal
  # rep("T1", 4),#9##Miller_Fisher
  rep("CTR", 3),#25##Control
  rep("CTR", 3),#26##Control
  rep("CTR", 3),#27##Control
  rep("CTR", 3),#28##Control
  rep("CTR", 3),#29##Control
  rep("T3", 3),#30##Unknown##1200
  rep("T3", 3),#31##Unknown##1200
  rep("T3", 3),#32##Unknown##1200
  rep("T3", 3),#33##Unknown##1200
  rep("T3", 3),#34##Unknown##1200
  rep("T3", 3),#35##Unknown##1200
  rep("T3", 3),#36##Unknown##1200
  rep("T3", 3),#37##Unknown##1200
  rep("T2", 3),#38##Axonal##60
  rep("T1", 3),#39##Axonal
  rep("T1", 3),#40##Demyelinating
  rep("T3", 3),#41##Demyelinating##1200
  rep("T1", 3),#42##Axonal
  rep("T1", 3),#43##Axonal
  rep("T2", 3),#44##Axonal##90
  rep("T2", 3),#45##Demyelinating##60
  rep("T1", 3),#46##Demyelinating
  rep("T2", 3),#47##Axonal##60
  rep("T1", 3)#48##Miller_Fisher
)


rownamescts <- cts[,1]
cts <- cts[-1]

#rownamescts <- gsub("\\..*$", "", rownames(cts))
rownames(cts) <- rownamescts
colnamescts <- colnames(cts)

cts <- cts[,-c(1:56)]

coldata <-  as.data.frame(cond)

colnames(cts) <- c("25_2", "25_3", "25_8", "26_2", "26_3", "26_8", "27_2", "27_3", "27_8", "28_2", "28_3", "28_8", "29_2", "29_3", "29_8", "30_2", "30_3", "30_8", "31_2", "31_3", "31_8",
                   "32_2", "32_3", "32_8", "33_2", "33_3", "33_8", "34_2", "34_3", "34_8", "35_2", "35_3", "35_8", "36_2", "36_3", "36_8", "37_2", "37_3", "37_8", "38_2", "38_3", "38_8", "39_2", "39_3", "39_8",
                   "40_2", "40_3", "40_8",  "41_2", "41_3", "41_8",  "42_2", "42_3", "42_8",  "43_2", "43_3", "43_8",
                   "44_2", "44_3", "44_8",  "45_2", "45_3", "45_8",  "46_2", "46_3", "46_8",  "47_2", "47_3", "47_8",
                   "48_2", "48_3", "48_8" )

#write.table(cts, file = "counts_GBS.csv",row.names=T, col.names = T,  sep=",")
#write.csv(cts, file = "counts_GBS.csv",row.names=T)
#teste <- read.table("counts_GBS.csv", header = T, stringsAsFactors = F, sep = ",")
rownames(coldata) <- colnames(cts)
colnames(coldata) <- c("condition")


type = c(#rep("Demyelinating", 4),#10
  #rep("Demyelinating", 4),#11
  # rep("Axonal", 4),#12
  # rep("Axonal", 4),#13
  # #rep("Axonal", 4),#14
  # rep("Demyelinating", 4),#15
  # #rep("Miller_Fisher", 4),#16
  # rep("Axonal", 4),#17
  # rep("Axonal", 4),#18
  # #rep("Axonal", 4),#19
  # #rep("Demyelinating", 4),#1
  # rep("Demyelinating", 4),#20
  # rep("Demyelinating", 4),#21
  # rep("Demyelinating", 4),#22
  # rep("Demyelinating", 4),#23
  # rep("Demyelinating", 4),#24
  # rep("Demyelinating", 4),#2
  # rep("Demyelinating", 4),#3
  # rep("Demyelinating", 4),#4
  # #rep("Axonal", 4),#5
  # #rep("Axonal", 4),#6
  # # rep("Axonal", 4),#7
  # #rep("Axonal", 4),#8
  # rep("Axonal", 4),#9
  rep("Control", 3),#25
  rep("Control", 3),#26
  rep("Control", 3),#27
  rep("Control", 3),#28
  rep("Control", 3),#29
  rep("Axonal", 3),#30
  rep("Axonal", 3),#31
  rep("Axonal", 3),#32
  rep("Demyelinating", 3),#33
  rep("Axonal", 3),#34
  rep("Axonal", 3),#35
  rep("Axonal", 3),#36
  rep("Demyelinating", 3),#37
  rep("Axonal", 3),#38
  rep("Axonal", 3),#39
  rep("Demyelinating", 3),#40
  rep("Demyelinating", 3),#41
  rep("Axonal", 3),#42
  rep("Axonal", 3),#43
  rep("Axonal", 3),#44
  rep("Demyelinating", 3),#45
  rep("Demyelinating", 3),#46
  rep("Axonal", 3),#47
  rep("Axonal", 3)#48
)


coldata <- cbind(coldata, type)
# coldata <- cbind(coldata, type, severity, days_to_walk)

all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds
# dds$samples <- factor(c(rep("10", 4),#10
#                         rep("11", 4),#11
#                         rep("12", 4),#12
#                         rep("13", 4),#13
#                         rep("14", 4),#14
#                         rep("15", 4),#15
#                         rep("16", 4),#16
#                         rep("17", 4),#17
#                         rep("18", 4),#18
#                         rep("19", 4),#19
#                         rep("1", 4),#1
#                         rep("20", 4),#20
#                         rep("21", 4),#21
#                         rep("22", 4),#22
#                         rep("23", 4),#23
#                         rep("24", 4),#24
#                         rep("2", 4),#2
#                         rep("3", 4),#3
#                         rep("4", 4),#4
#                         rep("5", 4),#5
#                         rep("6", 4),#6
#                         rep("7", 4),#7
#                         rep("8", 4),#8
#                         rep("9", 4),#9
#                         rep("25", 3),#25
#                         rep("26", 3),#26
#                         rep("27", 3),#27
#                         rep("28", 3),#28
#                         rep("29", 3),#29
#                         rep("30", 3),#30
#                         rep("31", 3),#31
#                         rep("32", 3),#32
#                         rep("33", 3),#33
#                         rep("34", 3),#34
#                         rep("35", 3),#35
#                         rep("36", 3),#36
#                         rep("37", 3),#37
#                         rep("38", 3),#38
#                         rep("39", 3),#39
#                         rep("40", 3),#40
#                         rep("41", 3),#41
#                         rep("42", 3),#42
#                         rep("43", 3),#43
#                         rep("44", 3),#44
#                         rep("45", 3),#45
#                         rep("46", 3),#46
#                         rep("47", 3),#47
#                         rep("48", 3)#48
# ))

featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
#dds$condition <- factor(dds$condition, levels = c("REC60","REC90","REC120","REC180", "REC240", "REC1000","GBS", "CTR"))
dds$condition <- factor(dds$condition, levels = c("T1","T2","T3", "CTR"))
dds$condition <- relevel(dds$condition, ref = "CTR")
dds$condition <- droplevels(dds$condition)
dds <- DESeq(dds)
res_t1 <- results(dds, alpha = 0.001)
summary(res_t1)

res_t1 <- results(dds, name="condition_T1_vs_CTR")
res_t1 <- results(dds, contrast=c("condition","T1","CTR"), alpha = 0.001)
summary(res_t1)
nrow(res_t1)
res_t1$ensembl <- rownames(res_t1)
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
                  filters = "ensembl_gene_id",
                  values = res_t1$ensembl,
                  mart = ensembl)
idx <- match(res_t1$ensembl, genemap$ensembl_gene_id )
#res$entrez <- genemap$entrezgene[idx]
res_t1$entrez <- genemap$entrezgene_id[idx]
res_t1$hgnc_symbol <- genemap$hgnc_symbol[idx]
summary(res_t1)
#res <- as.data.frame(res)
res01_t1 <- subset(res_t1, padj < 0.001 & abs(log2FoldChange)> 1.5)
summary(res01_t1)
nrow(res01_t1)
res01_t1_df <- as.data.frame(res01_t1)
res01_t1_df[res01_t1_df == ""] <- NA
sum(is.na(res01_t1_df))
res_genes_t1 <- na.omit(res01_t1_df)
sum(is.na(res_genes_t1))
sum(duplicated(res_genes_t1$hgnc_symbol))
sum(duplicated(res_genes_t1$entrez))
geneList_t1 <- res_genes_t1[,"log2FoldChange"]
names(geneList_t1) <- res_genes_t1$hgnc_symbol
# geneList_t1 <- subset(geneList_t1, select= -c(2))
# colnames(geneList_t1) <- NULL
# geneList_t1 <- as.numeric(geneList_t1)
# rownames(geneList_t1) <- geneList_t1 
up_res_genes_t1 <- subset(res_genes_t1, res_genes_t1$log2FoldChange > 0)
down_res_genes_t1 <- subset(res_genes_t1, res_genes_t1$log2FoldChange < 0)
# codes <- as.data.frame(res[, c("hgnc_symbol", "entrez")])
# dup_res <- duplicated(res)
res_t2 <- results(dds, name="condition_T2_vs_CTR")
res_t2 <- results(dds, contrast=c("condition","T2","CTR"), alpha = 0.001)
summary(res_t2)
res_t2$ensembl <- rownames(res_t2)
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap_t2 <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
                     filters = "ensembl_gene_id",
                     values = res_t2$ensembl,
                     mart = ensembl)
idx_t2 <- match(res_t2$ensembl, genemap_t2$ensembl_gene_id )
#res$entrez <- genemap$entrezgene[idx]
res_t2$entrez <- genemap_t2$entrezgene_id[idx_t2]
res_t2$hgnc_symbol <- genemap_t2$hgnc_symbol[idx_t2]
res01_t2 <- subset(res_t2, padj < 0.001 & abs(log2FoldChange)> 1.5)
summary(res01_t2)
nrow(res01_t2)
res01_t2_df <- as.data.frame(res01_t2)
res01_t2_df[res01_t2_df == ""] <- NA
sum(is.na(res01_t1_df))
res_genes_t2 <- na.omit(res01_t2_df)
sum(is.na(res_genes_t2))
sum(duplicated(res_genes_t2$hgnc_symbol))
sum(duplicated(res_genes_t2$entrez))
nrow(res_genes_t2)
geneList_t2 <- res_genes_t2[,"log2FoldChange"]
names(geneList_t2) <- res_genes_t2$hgnc_symbol
up_res_genes_t2 <- subset(res_genes_t2, res_genes_t2$log2FoldChange > 0)
down_res_genes_t2 <- subset(res_genes_t2, res_genes_t2$log2FoldChange < 0)
# res01_t1 <- subset(res_t1, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_t1)
# res_60 <- results(dds, name="condition_T60_vs_CTR")
# res_60 <- results(dds, contrast=c("condition","T60","CTR"))
# res_60$ensembl <- rownames(res_60)

# library( "biomaRt" )
# ensembl_60 = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap_60 <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                   filters = "ensembl_gene_id",
#                   values = res_60$ensembl,
#                   mart = ensembl_60)
# idx_60 <- match(res_60$ensembl, genemap$ensembl_gene_id )
# #res$entrez <- genemap$entrezgene[idx]
# res_60$hgnc_symbol <- genemap_60$hgnc_symbol[idx]
# res01_60 <- subset(res_60, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_60)

res_t3 <- results(dds, name="condition_T3_vs_CTR")
res_t3 <- results(dds, contrast=c("condition","T3","CTR"), alpha = 0.001)
summary(res_t3)
res_t3$ensembl <- rownames(res_t3)
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap_t3 <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
                     filters = "ensembl_gene_id",
                     values = res_t3$ensembl,
                     mart = ensembl)
idx_t3 <- match(res_t3$ensembl, genemap_t3$ensembl_gene_id )
#res$entrez <- genemap$entrezgene[idx]
res_t3$entrez <- genemap_t3$entrezgene_id[idx_t3]
res_t3$hgnc_symbol <- genemap_t3$hgnc_symbol[idx_t3]
#res_t2 <- as.data.frame(res_t2)
res01_t3 <- subset(res_t3, padj < 0.001 & abs(log2FoldChange)> 1.5)
summary(res01_t3)
nrow(res01_t3)
res01_t3_df <- as.data.frame(res01_t3)
res01_t3_df[res01_t3_df == ""] <- NA
sum(is.na(res01_t3_df))
res_genes_t3 <- na.omit(res01_t3_df)
sum(is.na(res_genes_t3))
sum(duplicated(res_genes_t3$hgnc_symbol))
sum(duplicated(res_genes_t3$entrez))
nrow(res_genes_t3)
geneList_t3 <- res_genes_t3[,"log2FoldChange"]
names(geneList_t3) <- res_genes_t3$hgnc_symbol
up_res_genes_t3 <- subset(res_genes_t3, res_genes_t3$log2FoldChange > 0)
down_res_genes_t3 <- subset(res_genes_t3, res_genes_t3$log2FoldChange < 0)

library("ggplot2")
library("RColorBrewer")
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("condition", "type"), ntop= 16000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= condition, shape=type)) +  geom_text(aes(label=name),vjust=2) +
  geom_point(size=4) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

library("genefilter")
topVarGenes <- head(order(-rowVars(assay(rld))),1000)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("condition","type")])
pheatmap(mat, annotation_col=df)

pheatmap(mat, cluster_rows=T, show_rownames = F, show_colnames = F,
         cluster_cols=T, annotation_col=df, annotation_legend = T)

