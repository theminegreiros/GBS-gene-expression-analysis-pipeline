 setwd("/home/themi/IOWA_notebook/counts")
#setwd("D:/IOWA_notebook/counts")
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

cond = c(rep("T2", 4),#10##Demyelinating
         rep("T1", 4),#11##Demyelinating
         rep("T2", 4),#12##Miller_Fisher
         rep("T1", 4),#13##Axonal
         rep("T2", 4),#14##Axonal
         rep("T1", 4),#15##Demyelinating
         rep("T1", 4),#16##Miller_Fisher
         rep("T2", 4),#17##Miller_Fisher
         rep("T1", 4),#18##Axonal
         rep("T2", 4),#19##Axonal
         rep("T1", 4),#1##Demyelinating
         rep("T2", 4),#20##Demyelinating
         rep("T1", 4),#21##Demyelinating
         rep("T2", 4),#22##Demyelinating
         rep("T2", 4),#23##Demyelinating
         rep("T1", 4),#24##Demyelinating
         rep("T2", 4),#2##Demyelinating
         rep("T1", 4),#3##Demyelinating
         rep("T2", 4),#4##Demyelinating
         rep("T1", 4),#5##Axonal
         rep("T2", 4),#6##Axonal
         rep("T1", 4),#7##Axonal
         rep("T2", 4),#8##Axonal
         rep("T1", 4),#9##Miller_Fisher
         rep("CTR", 3),#25##Control
         rep("CTR", 3),#26##Control
         rep("CTR", 3),#27##Control
         rep("CTR", 3),#28##Control
         rep("CTR", 3),#29##Control
         rep("T3", 3),#30##Unknown
         rep("T3", 3),#31##Unknown
         rep("REC", 3),#32##Unknown
         rep("T3", 3),#33##Unknown
         rep("T3", 3),#34##Unknown
         rep("T3", 3),#35##Unknown
         rep("T3", 3),#36##Unknown
         rep("T3", 3),#37##Unknown
         rep("REC", 3),#38##Axonal
         rep("GBS", 3),#39##Axonal
         rep("GBS", 3),#40##Demyelinating
         rep("T3", 3),#41##Demyelinating
         rep("GBS", 3),#42##Axonal
         rep("GBS", 3),#43##Axonal
         rep("REC", 3),#44##Axonal
         rep("REC", 3),#45##Demyelinating
         rep("GBS", 3),#46##Demyelinating
         rep("REC", 3),#47##Axonal
         rep("GBS", 3)#48##Miller_Fisher
)

rownamescts <- cts[,1]
cts <- cts[-1]
rownames(cts) <- rownamescts
colnamescts <- colnames(cts)
coldata <-  as.data.frame(cond)

colnames(cts) <- c("10_5", "10_6", "10_7", "10_8", "11_5","11_6","11_7","11_8","12_5","12_6", "12_7", "12_8", "13_5","13_6", "13_7", "13_8",
                   "14_5", "14_6","14_7","14_8", "15_5","15_6", "15_7", "15_8", "16_5","16_6","16_7","16_8", "17_5","17_6", "17_7", "17_8", 
                   "18_5","18_6", "18_7", "18_8", "19_5","19_6","19_7","19_8", "1_5","1_6","1_7","1_8", "20_5","20_6", "20_7", "20_8",
                   "21_5","21_6", "21_7", "21_8", "22_5","22_6", "22_7", "22_8", "23_5","23_6", "23_7", "23_8",
                   "24_5","24_6", "24_7", "24_8", "2_5","2_6", "2_7", "2_8", "3_5","3_6", "3_7", "3_8", "4_5","4_6", "4_7", "4_8", "5_5",
                   "5_6","5_7","5_8", "6_5", "6_6", "6_7","6_8", "7_5", "7_6","7_7", "7_8", "8_5", "8_6","8_7","8_8",
                   "9_5","9_6", "9_7", "9_8", "25_2", "25_3", "25_8", "26_2", "26_3", "26_8", "27_2", "27_3", "27_8",
                   "28_2", "28_3", "28_8", "29_2", "29_3", "29_8", "30_2", "30_3", "30_8", "31_2", "31_3", "31_8",
                   "32_2", "32_3", "32_8", "33_2", "33_3", "33_8", "34_2", "34_3", "34_8", "35_2", "35_3", "35_8",
                   "36_2", "36_3", "36_8", "37_2", "37_3", "37_8", "38_2", "38_3", "38_8", "39_2", "39_3", "39_8",
                   "40_2", "40_3", "40_8",  "41_2", "41_3", "41_8",  "42_2", "42_3", "42_8",  "43_2", "43_3", "43_8",
                   "44_2", "44_3", "44_8",  "45_2", "45_3", "45_8",  "46_2", "46_3", "46_8",  "47_2", "47_3", "47_8",
                   "48_2", "48_3", "48_8")

rownames(coldata) <- colnames(cts)
colnames(coldata) <- c("condition")

type = c(rep("Demyelinating", 4),#10
         rep("Demyelinating", 4),#11
         rep("Axonal", 4),#12
         rep("Axonal", 4),#13
         rep("Axonal", 4),#14
         rep("Demyelinating", 4),#15
         rep("Axonal", 4),#16
         rep("Axonal", 4),#17
         rep("Axonal", 4),#18
         rep("Axonal", 4),#19
         rep("Demyelinating", 4),#1
         rep("Demyelinating", 4),#20
         rep("Demyelinating", 4),#21
         rep("Demyelinating", 4),#22
         rep("Demyelinating", 4),#23
         rep("Demyelinating", 4),#24
         rep("Demyelinating", 4),#2
         rep("Demyelinating", 4),#3
         rep("Demyelinating", 4),#4
         rep("Axonal", 4),#5
         rep("Axonal", 4),#6
         rep("Axonal", 4),#7
         rep("Axonal", 4),#8
         rep("Axonal", 4),#9
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
coldata <- coldata[which(coldata$condition == "T1"|coldata$condition =="T2" |coldata$condition =="T3" | coldata$condition == "CTR"), ]
coldata <- coldata[, c(1,2)]
idx_cts <- rownames(coldata)

cts  <- cts[, idx_cts]
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds


featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
#dds$condition <- factor(dds$condition, levels = c("REC60","REC90","REC120","REC180", "REC240", "REC1000","GBS", "CTR"))
dds$condition <- factor(dds$condition, levels = c("T1","T2","T3", "CTR"))
dds$condition <- relevel(dds$condition, ref = "CTR")
dds$condition <- droplevels(dds$condition)
dds <- DESeq(dds)
res_t1 <- results(dds, alpha = 0.00001)
summary(res_t1)

library("ggplot2")
library("RColorBrewer")
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(rld, intgroup=c("condition", "type"), ntop= 16246, returnData=TRUE)
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
#res01_t1 <- subset(res_t1, padj < 0.001 & abs(log2FoldChange)> 1.5)
res01_t1 <- subset(res_t1, padj < 0.001)
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
res01_t2 <- subset(res_t2, padj < 0.001)
#res01_t2 <- subset(res_t2, padj < 0.001 & abs(log2FoldChange)> 1.5)
summary(res01_t2)
nrow(res01_t2)
res01_t2_df <- as.data.frame(res01_t2)
res01_t2_df[res01_t2_df == ""] <- NA
sum(is.na(res01_t2_df))
res_genes_t2 <- na.omit(res01_t2_df)
sum(is.na(res_genes_t2))
sum(duplicated(res_genes_t2$hgnc_symbol))
sum(duplicated(res_genes_t2$entrez))
nrow(res_genes_t2)
geneList_t2 <- res_genes_t2[,"log2FoldChange"]
names(geneList_t2) <- res_genes_t2$hgnc_symbol
up_res_genes_t2 <- subset(res_genes_t2, res_genes_t2$log2FoldChange > 0)
down_res_genes_t2 <- subset(res_genes_t2, res_genes_t2$log2FoldChange < 0)
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


library(org.Hs.eg.db)
library(clusterProfiler)
ego_t1 <- enrichGO(gene          = rownames(res01_t1),
                   keyType = "ENSEMBL",
                   universe   = rownames(cts),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pvalueCutoff  = 0.001,
                   #qvalueCutoff  = 0.05,
                   readable      = T)
head(ego_t1)

View(ego_t1@result)
ego_res_t1 <- ego_t1@result
ego_res_t1 <- subset(ego_res_t1, ego_res_t1$p.adjust <0.05)
unique(ego_res_t1$geneID)

dotplot(ego_t1, showCategory = 40)
barplot(ego_t1, showCategory = 10)
#color_palette(C("green", "red"))
cnetplot(ego_t1, showCategory = 10 , foldChange = geneList_t1)
heatplot(ego_t1, showCategory = 40, foldChange = geneList_t1)
# library(DOSE)
# data(geneList)
# de <- names(geneList)[abs(geneList) > 2]
# 
# edo <- enrichDGN(de)
# edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
# cnetplot(edox, foldChange=geneList)
# cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
# heatplot(ego_t1, foldChange = res_t1$log2FoldChange)


kk_t1 <- enrichKEGG(gene         = res_genes_t1$entrez,
                    organism     = 'hsa',
                    pvalueCutoff  = 0.1
)
head(kk_t1)

View(kk_t1@result)
barplot(kk_t1, drop=FALSE, showCategory=40)
cnetplot(kk_t1, showCategory = 10, foldChange = geneList_t1)
kk_t1x <- setReadable(kk_t1, 'org.Hs.eg.db', 'ENTREZID')
heatplot(kk_t1x, showCategory = 20 , foldChange = geneList_t1)



ego_t2 <- enrichGO(gene          = rownames(res01_t2),
                   keyType = "ENSEMBL",
                   universe   = rownames(cts),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pvalueCutoff  = 0.05,
                   readable      = T)


head(ego_t2)

View(ego_t2@result)
barplot(ego_t2, showCategory = 40)
cnetplot(ego_t2, showCategory = 10 , foldChange = geneList_t2)
heatplot(ego_t2, showCategory = 20 , foldChange = geneList_t2)
# library(DOSE)
# ego_60 <- enrichGO(gene          = rownames(res01_60),
#                 keyType = "ENSEMBL",
#                 universe   = rownames(cts),
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "BP",
#                 qvalueCutoff  = 0.05,
# #                 readable      = TRUE)
# head(ego_60)
# 
# View(ego_60@result)
# #egores <- ego@result
# #egores <- subset(egores, egores$p.adjust <0.05)
kk_t2 <- enrichKEGG(gene         = res_genes_t2$entrez,
                    organism     = 'hsa',
                    pvalueCutoff  = 0.05
)
head(kk_t2)
View(kk_t2@result)
barplot(kk_t2, drop=FALSE, showCategory=30)# 
kk_t2x <- setReadable(kk_t2, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kk_t2, showCategory = 10, foldChange = geneList_t2$log2FoldChange)
heatplot(kk_t2x, showCategory = 20 , foldChange = geneList_t2)
# 
# barplot(ego_60, showCategory = 20)
ego_t3 <- enrichGO(gene          = rownames(res01_t3),
                   keyType = "ENSEMBL",
                   universe   = rownames(cts),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)
head(ego_t3)

View(ego_t3@result)
#egores <- ego@result
#egores <- subset(egores, egores$p.adjust <0.05)
barplot(ego_t3, showCategory = 10)
cnetplot(ego_t3, showCategory = 10 , foldChange = geneList_t3)
heatplot(ego_t3, showCategory = 20 , foldChange = geneList_t3)
# library(DOSE)
# 
kk_t3 <- enrichKEGG(gene         = res_genes_t3$entrez,
                    organism     = 'hsa',
                    pvalueCutoff  = 0.01
)
head(kk_t3)
View(kk_t3@result)
barplot(kk_t3, drop=FALSE, showCategory=10)
emapplot(kk_t3)
kk_t3x <- setReadable(kk_t3, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kk_t3, showCategory = 10, foldChange = geneList_t3$log2FoldChange)
heatplot(kk_t3x, showCategory = 20 , foldChange = geneList_t3)
