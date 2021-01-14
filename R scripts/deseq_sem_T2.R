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


#cts <- cts[,-c(2:9, 38, 39, 40, 41)]
# cond = c(rep("REC60", 4),#10##Demyelinating
#          rep("GBS", 4),#11##Demyelinating
#          rep("REC120", 4),#12##Miller_Fisher
#          rep("GBS", 4),#13##Axonal
#          rep("REC60", 4),#14##Axonal
#          rep("GBS", 4),#15##Demyelinating
#          rep("GBS", 4),#16##Miller_Fisher
#          rep("REC60", 4),#17##Miller_Fisher
#          rep("GBS", 4),#18##Axonal
#          rep("REC120", 4),#19##Axonal
#          rep("GBS", 4),#1##Demyelinating
#          rep("REC180", 4),#20##Demyelinating
#          rep("GBS", 4),#21##Demyelinating
#          rep("REC240", 4),#22##Demyelinating
#          rep("REC60", 4),#23##Demyelinating
#          rep("GBS", 4),#24##Demyelinating
#          rep("REC60", 4),#2##Demyelinating
#          rep("GBS", 4),#3##Demyelinating
#          rep("REC120", 4),#4##Demyelinating
#          rep("GBS", 4),#5##Axonal
#          rep("REC60", 4),#6##Axonal
#          rep("GBS", 4),#7##Axonal
#          rep("REC90", 4),#8##Axonal
#          rep("GBS", 4),#9##Miller_Fisher
#          rep("CTR", 3),#25##Control
#          rep("CTR", 3),#26##Control
#          rep("CTR", 3),#27##Control
#          rep("CTR", 3),#28##Control
#          rep("CTR", 3),#29##Control
#          rep("REC1000", 3),#30##Unknown
#          rep("REC1000", 3),#31##Unknown
#          rep("REC1000", 3),#32##Unknown
#          rep("REC1000", 3),#33##Unknown
#          rep("REC1000", 3),#34##Unknown
#          rep("REC1000", 3),#35##Unknown
#          rep("REC1000", 3),#36##Unknown
#          rep("REC1000", 3),#37##Unknown
#          rep("REC60", 3),#38##Axonal
#          rep("GBS", 3),#39##Axonal
#          rep("GBS", 3),#40##Demyelinating
#          rep("REC1000", 3),#41##Demyelinating
#          rep("GBS", 3),#42##Axonal
#          rep("GBS", 3),#43##Axonal
#          rep("REC90", 3),#44##Axonal
#          rep("REC60", 3),#45##Demyelinating
#          rep("GBS", 3),#46##Demyelinating
#          rep("REC60", 3),#47##Axonal
#          rep("GBS", 3)#48##Miller_Fisher
# )

# cond = c(rep("REC60", 4),#10##Demyelinating
#          rep("GBS", 4),#11##Demyelinating
#          rep("REC120", 4),#12##Miller_Fisher
#          rep("GBS", 4),#13##Axonal
#          rep("REC60", 4),#14##Axonal
#          rep("GBS", 4),#15##Demyelinating
#          rep("GBS", 4),#16##Miller_Fisher
#          rep("REC60", 4),#17##Miller_Fisher
#          rep("GBS", 4),#18##Axonal
#          rep("GBS", 4),#1##Demyelinating
#          rep("REC180", 4),#20##Demyelinating
#          rep("GBS", 4),#21##Demyelinating
#          rep("REC240", 4),#22##Demyelinating
#          rep("REC60", 4),#23##Demyelinating
#          rep("GBS", 4),#24##Demyelinating
#          rep("REC60", 4),#2##Demyelinating
#          rep("GBS", 4),#3##Demyelinating
#          rep("REC120", 4),#4##Demyelinating
#          rep("GBS", 4),#5##Axonal
#          rep("REC60", 4),#6##Axonal
#          rep("GBS", 4),#7##Axonal
#          rep("REC90", 4),#8##Axonal
#          rep("GBS", 4),#9##Miller_Fisher
#          rep("CTR", 3),#25##Control
#          rep("CTR", 3),#26##Control
#          rep("CTR", 3),#27##Control
#          rep("CTR", 3),#28##Control
#          rep("CTR", 3),#29##Control
#          rep("REC1000", 3),#30##Unknown
#          rep("REC1000", 3),#31##Unknown
#          rep("REC1000", 3),#32##Unknown
#          rep("REC1000", 3),#33##Unknown
#          rep("REC1000", 3),#34##Unknown
#          rep("REC1000", 3),#35##Unknown
#          rep("REC1000", 3),#36##Unknown
#          rep("REC1000", 3),#37##Unknown
#          rep("REC60", 3),#38##Axonal
#          rep("GBS", 3),#39##Axonal
#          rep("GBS", 3),#40##Demyelinating
#          rep("REC1000", 3),#41##Demyelinating
#          rep("GBS", 3),#42##Axonal
#          rep("GBS", 3),#43##Axonal
#          rep("REC90", 3),#44##Axonal
#          rep("REC60", 3),#45##Demyelinating
#          rep("GBS", 3),#46##Demyelinating
#          rep("REC60", 3),#47##Axonal
#          rep("GBS", 3)#48##Miller_Fisher
# )
cond = c(#rep("REC60", 4),#10##Demyelinating
  #rep("GBS", 4),#11##Demyelinating
  rep("T2", 4),#12##Miller_Fisher##120
  rep("T1", 4),#13##Axonal
  #rep("REC60", 4),#14##Axonal
  rep("T1", 4),#15##Demyelinating
  #rep("GBS", 4),#16##Miller_Fisher
  rep("T2", 4),#17##Miller_Fisher##60
  rep("T1", 4),#18##Axonal
  #rep("REC120", 4),#19##Axonal
  #rep("GBS", 4),#1##Demyelinating
  rep("T2", 4),#20##Demyelinating##180
  rep("T1", 4),#21##Demyelinating
  rep("T2", 4),#22##Demyelinating##240
  rep("T2", 4),#23##Demyelinating##60
  rep("T1", 4),#24##Demyelinating
  rep("T2", 4),#2##Demyelinating##60
  rep("T1", 4),#3##Demyelinating
  rep("T2", 4),#4##Demyelinating##120
  #rep("GBS", 4),#5##Axonal
  #rep("REC60", 4),#6##Axonal
  #rep("GBS", 4),#7##Axonal
  #rep("REC90", 4),#8##Axonal
  rep("T1", 4),#9##Miller_Fisher
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

# cond = c(#rep("REC60", 4),#10##Demyelinating
#   #rep("GBS", 4),#11##Demyelinating
#   rep("T2", 4),#12##Miller_Fisher##120
#   rep("GBS", 4),#13##Axonal
#   #rep("REC60", 4),#14##Axonal
#   rep("GBS", 4),#15##Demyelinating
#   #rep("GBS", 4),#16##Miller_Fisher
#   rep("T1", 4),#17##Miller_Fisher##60
#   rep("GBS", 4),#18##Axonal
#   #rep("REC120", 4),#19##Axonal
#   #rep("GBS", 4),#1##Demyelinating
#   rep("T2", 4),#20##Demyelinating##180
#   rep("GBS", 4),#21##Demyelinating
#   rep("T2", 4),#22##Demyelinating##240
#   rep("T1", 4),#23##Demyelinating##60
#   rep("GBS", 4),#24##Demyelinating
#   rep("T1", 4),#2##Demyelinating##60
#   rep("GBS", 4),#3##Demyelinating
#   rep("T2", 4),#4##Demyelinating##120
#   #rep("GBS", 4),#5##Axonal
#   #rep("REC60", 4),#6##Axonal
#   #rep("GBS", 4),#7##Axonal
#   #rep("REC90", 4),#8##Axonal
#   rep("GBS", 4),#9##Miller_Fisher
#   rep("CTR", 3),#25##Control
#   rep("CTR", 3),#26##Control
#   rep("CTR", 3),#27##Control
#   rep("CTR", 3),#28##Control
#   rep("CTR", 3),#29##Control
#   rep("T3", 3),#30##Unknown##1200
#   rep("T3", 3),#31##Unknown##1200
#   rep("T3", 3),#32##Unknown##1200
#   rep("T3", 3),#33##Unknown##1200
#   rep("T3", 3),#34##Unknown##1200
#   rep("T3", 3),#35##Unknown##1200
#   rep("T3", 3),#36##Unknown##1200
#   rep("T3", 3),#37##Unknown##1200
#   rep("T1", 3),#38##Axonal##60
#   rep("GBS", 3),#39##Axonal
#   rep("GBS", 3),#40##Demyelinating
#   rep("T3", 3),#41##Demyelinating##1200
#   rep("GBS", 3),#42##Axonal
#   rep("GBS", 3),#43##Axonal
#   rep("T2", 3),#44##Axonal##90
#   rep("T1", 3),#45##Demyelinating##60
#   rep("GBS", 3),#46##Demyelinating
#   rep("T1", 3),#47##Axonal##60
#   rep("GBS", 3)#48##Miller_Fisher
# )
rownamescts <- cts[,1]
cts <- cts[-1]
#rownamescts <- gsub("\\..*$", "", rownames(cts))
rownames(cts) <- rownamescts
colnamescts <- colnames(cts)

#colnames(cts)<-gsub("^.*BAM.", "" , colnames(cts))
#colnames(cts)<-gsub(".bam", "" , colnames(cts))
#colnames(cts)<-gsub("_cliped50", "" , colnames(cts))
#colnames(cts)<-gsub("00", "", colnames(cts))

#write.table(cts, "COUNTSGBSRAW", sep = "\t", row.names = TRUE, col.names = NA)
#write.csv(cts, "COUNTSGBSRAW.csv")
coldata <-  as.data.frame(cond)

colnames(cts) <- c("12_5","12_6", "12_7", "12_8", "13_5","13_6", "13_7", "13_8", "15_5","15_6", "15_7", "15_8",
                   "17_5","17_6", "17_7", "17_8", "18_5","18_6", "18_7", "18_8", "20_5","20_6", "20_7", "20_8",
                   "21_5","21_6", "21_7", "21_8", "22_5","22_6", "22_7", "22_8", "23_5","23_6", "23_7", "23_8",
                   "24_5","24_6", "24_7", "24_8", "2_5","2_6", "2_7", "2_8", "3_5","3_6", "3_7", "3_8", "4_5","4_6", "4_7", "4_8",
                   "9_5","9_6", "9_7", "9_8", "25_2", "25_3", "25_8", "26_2", "26_3", "26_8", "27_2", "27_3", "27_8",
                   "28_2", "28_3", "28_8", "29_2", "29_3", "29_8", "30_2", "30_3", "30_8", "31_2", "31_3", "31_8",
                   "32_2", "32_3", "32_8", "33_2", "33_3", "33_8", "34_2", "34_3", "34_8", "35_2", "35_3", "35_8",
                   "36_2", "36_3", "36_8", "37_2", "37_3", "37_8", "38_2", "38_3", "38_8", "39_2", "39_3", "39_8",
                   "40_2", "40_3", "40_8",  "41_2", "41_3", "41_8",  "42_2", "42_3", "42_8",  "43_2", "43_3", "43_8",
                   "44_2", "44_3", "44_8",  "45_2", "45_3", "45_8",  "46_2", "46_3", "46_8",  "47_2", "47_3", "47_8",
                   "48_2", "48_3", "48_8")

write.table(cts, file = "counts_GBS.csv",row.names=T, col.names = T,  sep=",")
write.csv(cts, file = "counts_GBS.csv",row.names=T)
teste <- read.table("counts_GBS.csv", header = T, stringsAsFactors = F, sep = ",")
rownames(coldata) <- colnames(cts)
colnames(coldata) <- c("condition")
# type = c(rep("Demyelinating", 4),#10
#          rep("Demyelinating", 4),#11
#          rep("Miller_Fisher", 4),#12
#          rep("Axonal", 4),#13
#          rep("Axonal", 4),#14
#          rep("Demyelinating", 4),#15
#          rep("Miller_Fisher", 4),#16
#          rep("Miller_Fisher", 4),#17
#          rep("Axonal", 4),#18
#          rep("Axonal", 4),#19
#          rep("Demyelinating", 4),#1
#          rep("Demyelinating", 4),#20
#          rep("Demyelinating", 4),#21
#          rep("Demyelinating", 4),#22
#          rep("Demyelinating", 4),#23
#          rep("Demyelinating", 4),#24
#          rep("Demyelinating", 4),#2
#          rep("Demyelinating", 4),#3
#          rep("Demyelinating", 4),#4
#          rep("Axonal", 4),#5
#          rep("Axonal", 4),#6
#          rep("Axonal", 4),#7
#          rep("Axonal", 4),#8
#          rep("Miller_Fisher", 4),#9
#          rep("Control", 3),#25
#          rep("Control", 3),#26
#          rep("Control", 3),#27
#          rep("Control", 3),#28
#          rep("Control", 3),#29
#          rep("Axonal", 3),#30
#          rep("Miller_Fisher", 3),#31
#          rep("Unknown", 3),#32
#          rep("Demyelinating", 3),#33
#          rep("Axonal", 3),#34
#          rep("Axonal", 3),#35
#          rep("Axonal", 3),#36
#          rep("Demyelinating", 3),#37
#          rep("Axonal", 3),#38
#          rep("Axonal", 3),#39
#          rep("Demyelinating", 3),#40
#          rep("Demyelinating", 3),#41
#          rep("Axonal", 3),#42
#          rep("Axonal", 3),#43
#          rep("Axonal", 3),#44
#          rep("Demyelinating", 3),#45
#          rep("Demyelinating", 3),#46
#          rep("Axonal", 3),#47
#          rep("Miller_Fisher", 3)#48
# )
# type = c(rep("Demyelinating", 4),#10
#          rep("Demyelinating", 4),#11
#          rep("Miller_Fisher", 4),#12
#          rep("Axonal", 4),#13
#          rep("Axonal", 4),#14
#          rep("Demyelinating", 4),#15
#          rep("Miller_Fisher", 4),#16
#          rep("Miller_Fisher", 4),#17
#          rep("Axonal", 4),#18
#          rep("Demyelinating", 4),#1
#          rep("Demyelinating", 4),#20
#          rep("Demyelinating", 4),#21
#          rep("Demyelinating", 4),#22
#          rep("Demyelinating", 4),#23
#          rep("Demyelinating", 4),#24
#          rep("Demyelinating", 4),#2
#          rep("Demyelinating", 4),#3
#          rep("Demyelinating", 4),#4
#          rep("Axonal", 4),#5
#          rep("Axonal", 4),#6
#          rep("Axonal", 4),#7
#          rep("Axonal", 4),#8
#          rep("Miller_Fisher", 4),#9
#          rep("Control", 3),#25
#          rep("Control", 3),#26
#          rep("Control", 3),#27
#          rep("Control", 3),#28
#          rep("Control", 3),#29
#          rep("Axonal", 3),#30
#          rep("Miller_Fisher", 3),#31
#          rep("Miller_Fisher", 3),#32
#          rep("Demyelinating", 3),#33
#          rep("Axonal", 3),#34
#          rep("Axonal", 3),#35
#          rep("Axonal", 3),#36
#          rep("Demyelinating", 3),#37
#          rep("Axonal", 3),#38
#          rep("Axonal", 3),#39
#          rep("Demyelinating", 3),#40
#          rep("Demyelinating", 3),#41
#          rep("Axonal", 3),#42
#          rep("Axonal", 3),#43
#          rep("Axonal", 3),#44
#          rep("Demyelinating", 3),#45
#          rep("Demyelinating", 3),#46
#          rep("Axonal", 3),#47
#          rep("Miller_Fisher", 3)#48
# )

type = c(#rep("Demyelinating", 4),#10
  #rep("Demyelinating", 4),#11
  rep("Axonal", 4),#12
  rep("Axonal", 4),#13
  #rep("Axonal", 4),#14
  rep("Demyelinating", 4),#15
  #rep("Miller_Fisher", 4),#16
  rep("Axonal", 4),#17
  rep("Axonal", 4),#18
  #rep("Axonal", 4),#19
  #rep("Demyelinating", 4),#1
  rep("Demyelinating", 4),#20
  rep("Demyelinating", 4),#21
  rep("Demyelinating", 4),#22
  rep("Demyelinating", 4),#23
  rep("Demyelinating", 4),#24
  rep("Demyelinating", 4),#2
  rep("Demyelinating", 4),#3
  rep("Demyelinating", 4),#4
  #rep("Axonal", 4),#5
  #rep("Axonal", 4),#6
  # rep("Axonal", 4),#7
  #rep("Axonal", 4),#8
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
# severity = c(#rep("0", 4),#10##Demyelinating
#          #rep("0", 4),#11##Demyelinating
#          rep("0", 4),#12##Miller_Fisher
#          rep("1", 4),#13##Axonal
#          #rep("1", 4),#14##Axonal
#          rep("1", 4),#15##Demyelinating
#          #rep("1", 4),#16##Miller_Fisher
#          rep("1", 4),#17##Miller_Fisher
#          rep("1", 4),#18##Axonal
#          #rep("1", 4),#1##Demyelinating
#          rep("1", 4),#20##Demyelinating
#          rep("1", 4),#21##Demyelinating
#          rep("1", 4),#22##Demyelinating
#          rep("0", 4),#23##Demyelinating
#          rep("0", 4),#24##Demyelinating
#          rep("0", 4),#2##Demyelinating
#          rep("1", 4),#3##Demyelinating
#          rep("1", 4),#4##Demyelinating
#          #rep("1", 4),#5##Axonal
#          #rep("1", 4),#6##Axonal
#          #rep("1", 4),#7##Axonal
#          #rep("1", 4),#8##Axonal
#          rep("0", 4),#9##Miller_Fisher
#          rep("CTR", 3),#25##Control
#          rep("CTR", 3),#26##Control
#          rep("CTR", 3),#27##Control
#          rep("CTR", 3),#28##Control
#          rep("CTR", 3),#29##Control
#          rep("1", 3),#30##Unknown
#          rep("1", 3),#31##Unknown
#          rep("1", 3),#32##Unknown
#          rep("0", 3),#33##Unknown
#          rep("1", 3),#34##Unknown
#          rep("1", 3),#35##Unknown
#          rep("1", 3),#36##Unknown
#          rep("0", 3),#37##Unknown
#          rep("1", 3),#38##Axonal
#          rep("1", 3),#39##Axonal
#          rep("0", 3),#40##Demyelinating
#          rep("0", 3),#41##Demyelinating
#          rep("1", 3),#42##Axonal
#          rep("1", 3),#43##Axonal
#          rep("1", 3),#44##Axonal
#          rep("0", 3),#45##Demyelinating
#          rep("0", 3),#46##Demyelinating
#          rep("1", 3),#47##Axonal
#          rep("1", 3)#48##Miller_Fisher
# )
# 
# days_to_walk  = c(#rep("0", 4),#10##Demyelinating
#   #rep("0", 4),#11##Demyelinating
#   rep("0", 4),#12##Miller_Fisher
#   rep("0", 4),#13##Axonal
#   #rep("1", 4),#14##Axonal
#   rep("103", 4),#15##Demyelinating_103
# 
#   #rep("1", 4),#16##Miller_Fisher
#   rep("10", 4),#17##Miller_Fisher_10
#   rep("546", 4),#18##Axonal_546
#   #rep("1", 4),#1##Demyelinating
#   rep("103", 4),#20##Demyelinating
#   rep("25", 4),#21##Demyelinating
#   rep("25", 4),#22##Demyelinating
#   rep("0", 4),#23##Demyelinating
#   rep("0", 4),#24##Demyelinating
#   rep("0", 4),#2##Demyelinating
#   rep("60", 4),#3##Demyelinating
#   rep("60", 4),#4##Demyelinating
#   #rep("1", 4),#5##Axonal
#   #rep("1", 4),#6##Axonal
#   #rep("1", 4),#7##Axonal
#   #rep("1", 4),#8##Axonal
#   rep("0", 4),#9##Miller_Fisher
#   rep("CTR", 3),#25##Control
#   rep("CTR", 3),#26##Control
#   rep("CTR", 3),#27##Control
#   rep("CTR", 3),#28##Control
#   rep("CTR", 3),#29##Control
#   rep("100", 3),#30##Unknown
#   rep("103", 3),#31##Unknown
#   rep("0", 3),#32##Unknown
#   rep("0", 3),#33##Unknown
#   rep("0", 3),#34##Unknown
#   rep("45", 3),#35##Unknown
#   rep("546", 3),#36##Unknown
#   rep("0", 3),#37##Unknown
#   rep("100", 3),#38##Axonal
#   rep("546", 3),#39##Axonal
#   rep("0", 3),#40##Demyelinating
#   rep("0", 3),#41##Demyelinating
#   rep("100", 3),#42##Axonal
#   rep("45", 3),#43##Axonal
#   rep("45", 3),#44##Axonal
#   rep("0", 3),#45##Demyelinating
#   rep("0", 3),#46##Demyelinating
#   rep("0", 3),#47##Axonal
#   rep("10", 3)#48##Miller_Fisher
# )

coldata <- cbind(coldata, type)
# coldata <- cbind(coldata, type, severity, days_to_walk)
coldata <- coldata[which(coldata$condition == "T1"|coldata$condition =="T3" | coldata$condition == "CTR"), ]
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
keep <- rowSums(counts(dds)) >= 20
dds <- dds[keep,]
#dds$condition <- factor(dds$condition, levels = c("REC60","REC90","REC120","REC180", "REC240", "REC1000","GBS", "CTR"))
dds$condition <- factor(dds$condition, levels = c("T1","T3", "CTR"))
dds$condition <- relevel(dds$condition, ref = "CTR")
dds$condition <- droplevels(dds$condition)
dds <- DESeq(dds)
res_t1 <- results(dds, alpha = 0.001)
summary(res_t1)
# 
# dds$condition <- relevel(dds$condition, ref = "GBS")
# dds$condition <- droplevels(dds$condition)
# dds <- DESeq(dds)
#res
# #normalizando os counts para usar no cemitool
# dds <- estimateSizeFactors(dds)
# sizeFactors(dds)
# normalized_counts <- counts(dds, normalized=TRUE)
# normalized_counts <- as.data.frame(normalized_counts)
# normalized_counts$ensembl <- rownames(normalized_counts)
# library( "biomaRt" )
# ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                   filters = "ensembl_gene_id",
#                   values = normalized_counts$ensembl,
#                   mart = ensembl)
# 
# idx <- match(normalized_counts$ensembl, genemap$ensembl_gene_id )
# #res$entrez <- genemap$entrezgene[idx]
# normalized_counts$hgnc_symbol <- genemap$hgnc_symbol[idx]
# normalized_counts[normalized_counts==""] <- NA
# sum(is.na(normalized_counts))
# normalized_counts <- na.omit(normalized_counts) 
# sum(is.na(normalized_counts))
# codes <- normalized_counts[, c("ensembl", "hgnc_symbol")] 
# rownames(normalized_counts)  <- normalized_counts$hgnc_symbol
# falses <- duplicated(normalized_counts$ensembl)
# falses_2 <- duplicated(normalized_counts$hgnc_symbol)
# remove <- c("ENSG00000284770")
# normalized_counts_2 <- normalized_counts[!(rownames(normalized_counts) %in% remove),]
# rownames(normalized_counts_2) <- normalized_counts_2$hgnc_symbol 
# cols.dont.want <- c("ensembl","hgnc_symbol")
# normalized_counts_2 <- normalized_counts_2[, ! names(normalized_counts_2) %in% cols.dont.want, drop = F] 
#   #normalized_counts_2[,-c(209, 210)]
# #cemitool
# library("CEMiTool")
# cem <- cemitool(normalized_counts_2, set_beta= 14)
# cem <- diagnostic_report(cem, force=T)
# get_cemitool_r2_beta(cem, eps = 0.1)
# cem <- cemitool(expr=normalized_counts_2)
# generate_report(cem, force = T)
# nmodules(cem)
# head(module_genes(cem))
# generate_report(cem, force=T)
# write_files(cem, force = T)
# save_plots(cem, "all", force =T)
# sample_annot <- coldata
# sample_annot$SampleName <-rownames(sample_annot)
# sample_annot <- sample_annot[,c("SampleName", "condition")]
# colnames(sample_annot) <- c("SampleName", "Class")
# cem <- cemitool(normalized_counts_2, sample_annot, set_beta= 14)
# sample_annotation(cem, 
#                   sample_name_column="SampleName", 
#                   class_column="Class") <- sample_annot
# # generate heatmap of gene set enrichment analysis
# cem <- mod_gsea(cem)
# cem <- plot_gsea(cem)
# show_plot(cem, "gsea")
# # plot gene expression within each module
# cem <- plot_profile(cem)
# plots <- show_plot(cem, "profile")
# plots[2]
# # Run CEMiTool applying Variance Stabilizing Transformation to data
# cem <- cemitool(expr=normalized_counts, apply_vst=TRUE)
# # Run CEMiTool with additional processing messages
# cem <- cemitool(expr=expr0, verbose=TRUE)

# normalized_counts_2 <- subset(normalized_counts_2, select = -c("ensembl","hgnc_symbol"))

#expressão diferencial
library("ggplot2")
library("RColorBrewer")
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), ntop= 10000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color= condition, shape=type)) +  geom_text(aes(label=name),vjust=2) +
  geom_point(size=4) + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
# library(plotly)
# plot_ly(pcaData) 
# # +  geom_text(aes(label=name),vjust=2) +
#   geom_point(size=4) + 
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   coord_fixed()






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
# up_res_genes_t1 <- subset(res_genes_t1, res_genes_t1$log2FoldChange > 0)
# down_res_genes_t1 <- subset(res_genes_t1, res_genes_t1$log2FoldChange < 0)
# # codes <- as.data.frame(res[, c("hgnc_symbol", "entrez")])
# # dup_res <- duplicated(res)
# res_t2 <- results(dds, name="condition_T2_vs_CTR")
# res_t2 <- results(dds, contrast=c("condition","T2","CTR"), alpha = 0.001)
# summary(res_t2)
# res_t2$ensembl <- rownames(res_t2)
# library( "biomaRt" )
# ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap_t2 <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id"),
#                      filters = "ensembl_gene_id",
#                      values = res_t2$ensembl,
#                      mart = ensembl)
# idx_t2 <- match(res_t2$ensembl, genemap_t2$ensembl_gene_id )
# #res$entrez <- genemap$entrezgene[idx]
# res_t2$entrez <- genemap_t2$entrezgene_id[idx_t2]
# res_t2$hgnc_symbol <- genemap_t2$hgnc_symbol[idx_t2]
# res01_t2 <- subset(res_t2, padj < 0.001 & abs(log2FoldChange)> 1.5)
# summary(res01_t2)
# nrow(res01_t2)
# res01_t2_df <- as.data.frame(res01_t2)
# res01_t2_df[res01_t2_df == ""] <- NA
# sum(is.na(res01_t1_df))
# res_genes_t2 <- na.omit(res01_t2_df)
# sum(is.na(res_genes_t2))
# sum(duplicated(res_genes_t2$hgnc_symbol))
# sum(duplicated(res_genes_t2$entrez))
# nrow(res_genes_t2)
# geneList_t2 <- res_genes_t2[,"log2FoldChange"]
# names(geneList_t2) <- res_genes_t2$hgnc_symbol
# up_res_genes_t2 <- subset(res_genes_t2, res_genes_t2$log2FoldChange > 0)
# down_res_genes_t2 <- subset(res_genes_t2, res_genes_t2$log2FoldChange < 0)
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
# downres02 <- subset(res02, res02$log2FoldChange < 0)
library("genefilter")
topVarGenes <- head(order(-rowVars(assay(rld))),1000)
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("condition","type")])
pheatmap(mat, annotation_col=df)

pheatmap(mat, cluster_rows=T, show_rownames = T, show_colnames = F,
         cluster_cols=T, annotation_col=df, annotation_legend = T)



# res01_t2 <- subset(res_t2, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_t2)

# library( "biomaRt" )
# ensembl_1200 = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap_1200 <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                      filters = "ensembl_gene_id",
#                      values = res_1200$ensembl,
#                      mart = ensembl_1200)
# idx_60 <- match(res_1200$ensembl, genemap$ensembl_gene_id )
#res$entrez <- genemap$entrezgene[idx]
# res_60$hgnc_symbol <- genemap_60$hgnc_symbol[idx]
# res01_60 <- subset(res_60, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_60)
# 
# res_60_gbs <- results(dds, name="condition_T60_vs_GBS")
# res_60_gbs <- results(dds, contrast=c("condition","T60","GBS"))
# nrow(res_60_gbs)
# res01_60_gbs <- subset(res_60_gbs, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_60_gbs)
# 
# res_1200_gbs <- results(dds, name="condition_T1200_vs_GBS")
# res_1200_gbs <- results(dds, contrast=c("condition","T1200","GBS"))
# nrow(res_1200_gbs)
# res01_1200_gbs <- subset(res_1200_gbs, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_1200_gbs)
# res_120_gbs <- results(dds, name="condition_T120_vs_GBS")
# res_120_gbs <- results(dds, contrast=c("condition","T120","GBS"))
# nrow(res_120_gbs)
# 
# res01_120_gbs <- subset(res_120_gbs, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_120_gbs)
# 
# res_240_gbs <- results(dds, name="condition_T240_vs_GBS")
# res_240_gbs <- results(dds, contrast=c("condition","T240","GBS"))
# nrow(res_240_gbs)
# 
# res01_240_gbs <- subset(res_240_gbs, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_240_gbs)
# 
# class(res01_240_gbs)
# df_120_gbs <- as.data.frame(res01_120_gbs)
# df_120_gbs$ID <- rownames(df_120_gbs)
# df_120_gbs <- df_120_gbs[,c(2, 7)]
# df_240_gbs <- as.data.frame(res01_240_gbs)
# df_240_gbs$ID <- rownames(df_240_gbs)
# df_240_gbs <- df_240_gbs[,c(2, 7)]
# head( res[ order( res$log2FoldChange ), ] )
# plotMA(res, ylim = c(-4, 4) )
# #res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
# library( "biomaRt" )
# ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
#                   filters = "ensembl_gene_id",
#                   values = res$ensembl,
#                   mart = ensembl )
# idx <- match( res$ensembl, genemap$ensembl_gene_id )
# res$entrez <- genemap$entrezgene[ idx ]
# res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]

# 
# res_90 <- results(dds, name="condition_REC90_vs_CTR")
# res_90 <- results(dds, contrast=c("condition","REC90","CTR"))
# res_120 <- results(dds, name="condition_T120_vs_CTR")
# res_120 <- results(dds, contrast=c("condition","T120","CTR"))
# res01_120 <- subset(res_120, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_120)
# res_240 <- results(dds, name="condition_T240_vs_CTR")
# res_240 <- results(dds, contrast=c("condition","T240","CTR"))
# res01_240 <- subset(res_240, padj < 0.001 & abs(log2FoldChange)> 2)
# nrow(res01_240)
# res_1200 <- results(dds, name="condition_T1200_vs_CTR")
# res_1200 <- results(dds, contrast=c("condition","T1200","CTR"))
# res01_1200 <- subset(res_1200, padj < 0.001 & abs(log2FoldChange)> 1.5)
# nrow(res01_1200)
# # head(res,4)
#  res01 <- subset(res, padj < 0.001 & abs(log2FoldChange)> 1.5)
#  #res01_60 <- subset(res_60, padj < 0.001 & abs(log2FoldChange)> 1.5)
# #nrow(res01_60)
# 
# 
# res01_240 <- subset(res_240, padj < 0.001 & abs(log2FoldChange)> 2)
# nrow(res01_240)
# 
# 
#  l# res01 <- subset(res, padj < 0.001)
# # plotMA( res01, ylim = c(-4, 4) )
#  nrow(res01)
# head( res01[ order( res01$log2FoldChange ), ] )
#plotDispEsts( dds, ylim = c(1e-6, 1e3) )
#res01 <- as.data.frame(res01)
# upres01 <- subset(res01, res01$log2FoldChange > 0)
# downres01 <- subset(res01, res01$log2FoldChange < 0)
# plotMA( upres01, ylim = c(-4, 4) )

# hist( res$pvalue, breaks=20, col="grey" )
# write.csv( as.data.frame(res01), header= T,  file="res01.csv" )
# write.table(resdata1, "resdata1", sep = "\t", row.names = TRUE, col.names = NA)
# resdata1 <- as.data.frame(res01)
# resdata <- as.data.frame(res)
# cdata1 <- rownames(res01)
# rownames(resdata1) <- NULL
# resdata1 <- cbind(cdata1,resdata1)
# resdata1 <- resdata1[,c(1, 3, 7)]
# all(rownames(res01) %in% resdata1$cdata1)
# all(rownames(res01) == resdata1$cdata1)
# write.table(resdata1, "resdata1", sep = "\t", row.names = F, col.names = F)
#
# #resdata[,1]
# #resdata <- resdata[,c(2, 6)]
# write.table(resdata, "resdata", sep = "\t", row.names = TRUE, col.names = NA)
#
# View(resdata1)
# resdata1[,1]
# idxres <- resdata1[,c(2, 6)]
# resdata2 <- idxres
# write.table(resdata2, "resdata2", sep = "\t", row.names = TRUE, col.names = NA)
# library( "biomaRt" )
# ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap <- getBM( attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                   filters = "ensembl_gene_id",
#                   values = res$ensembl,
#                   mart = ensembl)
# idx <- match(res$ensembl, genemap$ensembl_gene_id )
# #res$entrez <- genemap$entrezgene[idx]
# res$hgnc_symbol <- genemap$hgnc_symbol[idx]
# #res02 <- subset(res2, padj < 0.05 & abs(log2FoldChange)>1)
# nrow(res02) 

# library("ggplot2")
# library("RColorBrewer")
# vsd <- vst(dds, blind=FALSE)
# rld <- rlog(dds, blind=FALSE)
# # pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), ntop= 21000, returnData=TRUE)
# # percentVar <- round(100 * attr(pcaData, "percentVar"))
# # ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) + 
# #   geom_point(size=4) + scale_color_manual(values = c("indianred1","orange1","limegreen", "gray19", "aquamarine4","red4")) +
# #   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
# #   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
# #   coord_fixed()
# 
# pcaData <- plotPCA(rld, intgroup=c("condition", "type"), ntop= 16000, returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color= condition, shape=type)) + 
#   geom_point(size=4) + 
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
#   coord_fixed()
# scale_color_manual(values = c("indianred1","orange1","limegreen", "gray19", "aquamarine4","red4"))
# vsd <- vst(dds, blind=FALSE)
# pcaData <- plotPCA(vsd, intgroup=c("condition", "type"), ntop= 21000, returnData=TRUE)
# percentVar <- round(100 * attr(pcaData, "percentVar"))
# ggplot(pcaData, aes(PC1, PC2, color=condition, shape=type)) +
#   geom_point(size=2) +
#   xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black", nudge_x = 2, nudge_y = -2, size = 1) +
#   coord_fixed()
# library("limma")
# #VEEN DIAGRAM
# res_GBSvsCTR <- res01
# res_GBSvsCTR.genes <- row.names(res_GBSvsCTR)
# res_T60vsCTR <- res01_60
# res_T60vsCTR.genes <-row.names(res_T60vsCTR)
# res_T120vsCTR <- res01_120
# res_T120vsCTR.genes <- row.names(res_T120vsCTR)
# res_T240vsCTR <- res01_240
# res_T240vsCTR.genes <- row.names(res_T240vsCTR)
# res_T1200vsCTR <- res01_1200
# res_T1200vsCTR.genes <-row.names(res_T1200vsCTR)
# res_T240vsGBS <- res01_240_gbs
# res_T240vsGBS.genes <- row.names(res_T240vsGBS)
# res_T120vsGBS <- res_120_gbs
# res_T120vsGBS.genes <- row.names(res_T120vsGBS) 
# res_T1200vsGBS <- res_1200_gbs
# res_T1200vsGBS.genes <- row.names(res_T1200vsGBS)
# # Combining the two above..
# comb <- c(res_GBSvsCTR.genes, res_T60vsCTR.genes,res_T120vsCTR.genes, res_T240vsCTR.genes, res_T1200vsCTR.genes,res_T120vsGBS.genes, res_T240vsGBS.genes, res_T1200vsGBS.genes)
# 
# # Comparing comb with the above two
# res_GBSvsCTR.genes.2 <- comb %in% res_GBSvsCTR.genes
# res_T120vsCTR.genes.2 <- comb %in% res_T120vsCTR.genes 
# res_T240vsCTR.genes.2 <- comb %in% res_T240vsCTR.genes
# res_T1200vsCTR.genes.2 <- comb %in% res_T1200vsCTR.genes 
# res_T120vsGBS.genes.2 <- comb %in% res_T120vsGBS.genes
# res_T240vsGBS.genes.2 <- comb %in% res_T240vsGBS.genes 
# res_T1200vsGBS.genes.2 <- comb %in% res_T1200vsGBS.genes

# # Generating venn counts to plot venn diagram
# counts_33 <- cbind(res_SGBvsCTR.genes.2, res_T120vsCTR.genes.2, res_T240vsCTR.genes.2, res_T1200vsCTR.genes.2, res_T120vsGBS.genes.2, res_T240vsGBS.genes.2, res_T1200vsGBS.genes.2 )
# results.33 <- vennCounts(counts_33)
# vennDiagram(results.33, cex = 1,names = c("SGBvsCTR","T120vsCTR", "T240vsCTR","T1200vsCTR", "T120vSGB", "T1200vsSGB"), circle.col = c("red", "blue","yellow","green", "purple", "grey" ))
# 
# counts_43 <- cbind(res_GBSvsCTR.genes.2, res_T120vsCTR.genes.2, res_T1200vsCTR.genes.2, res_T120vsGBS.genes.2, res_T1200vsGBS.genes.2 )
# results.43 <- vennCounts(counts_43)
# vennDiagram(results.43, cex = 1,names = c("SGBvsCTR","T120vsCTR","T1200vsCTR", "T120vSGB", "T1200vsSGB"), circle.col = c("red", "blue","yellow","green", "grey" ))



library(org.Hs.eg.db)
library(clusterProfiler)
ego_t1 <- enrichGO(gene          = rownames(res01_t1),
                   keyType = "ENSEMBL",
                   universe   = rownames(cts),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pvalueCutoff  = 0.05,
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
heatplot(ego_t1, showCategory = 10, foldChange = geneList_t1)
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
barplot(kk_t1, drop=FALSE, showCategory=10)
cnetplot(kk_t1, showCategory = 10, foldChange = geneList_t1)
kk_t1x <- setReadable(kk_t1, 'org.Hs.eg.db', 'ENTREZID')
heatplot(kk_t1x, showCategory = 20 , foldChange = geneList_t1)



ego_t2 <- enrichGO(gene          = rownames(res01_t2),
                   keyType = "ENSEMBL",
                   universe   = rownames(cts),
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pvalueCutoff  = 0.1,
                   readable      = T)


head(ego_t2)

View(ego_t2@result)
barplot(ego_t2, showCategory = 10)
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
barplot(kk_t2, drop=FALSE, showCategory=10)# 
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


plotCounts(dds, gene= "ENSG00000162772", intgroup= "condition", normalized = T)
data <- plotCounts(dds, gene="ENSG00000115590", intgroup="condition", returnData=TRUE, normalized = T)
ggplot(data, aes(x=condition, y= count, fill=condition)) +
  scale_y_log10() + 
  geom_dotplot(binaxis="y", stackdir="center")
# ego_90 <- enrichGO(gene          = rownames(res01_90),
#                    keyType = "ENSEMBL",
#                    universe   = rownames(cts),
#                    OrgDb         = org.Hs.eg.db,
#                    ont           = "BP",
#                    qvalueCutoff  = 0.05,
#                    readable      = TRUE)
# head(ego_90)
# 
# View(ego_90@result)
# #egores <- ego@result
# #egores <- subset(egores, egores$p.adjust <0.05)
# 
# 
# barplot(ego_90, showCategory = 20)

# ego_120 <- enrichGO(gene          = rownames(res01_120),
#                    keyType = "ENSEMBL",
#                    universe   = rownames(cts),
#                    OrgDb         = org.Hs.eg.db,
#                    ont           = "BP",
#                    qvalueCutoff  = 0.05,
#                    readable      = TRUE)
# head(ego_120)
# 
# View(ego_120@result)
# #egores <- ego@result
# #egores <- subset(egores, egores$p.adjust <0.05)
# 
# 
# barplot(ego_120, showCategory = 20)
# 
# ego_240 <- enrichGO(gene          = rownames(res01_240),
#                     keyType = "ENSEMBL",
#                     universe   = rownames(cts),
#                     OrgDb         = org.Hs.eg.db,
#                     ont           = "BP",
#                     qvalueCutoff  = 0.05,
#                     readable      = TRUE)
# head(ego_240)
# 
# View(ego_2400@result)
# #egores <- ego@result
# #egores <- subset(egores, egores$p.adjust <0.05)
# 
# 
# barplot(ego_240, showCategory = 20)
# 
# 
# ego_1200 <- enrichGO(gene          = rownames(res01_1200),
#                     keyType = "ENSEMBL",
#                     universe   = rownames(cts),
#                     OrgDb         = org.Hs.eg.db,
#                     ont           = "BP",
#                     qvalueCutoff  = 0.05,
#                     readable      = TRUE)
# head(ego_1200)
# 
# View(ego_1200@result)
# #egores <- ego@result
# #egores <- subset(egores, egores$p.adjust <0.05)
# barplot(ego_1200, showCategory = 20)

# ego_1200 <- enrichGO(gene          = rownames(res01_1200),
#                      keyType = "ENSEMBL",
#                      universe   = rownames(cts),
#                      OrgDb         = org.Hs.eg.db,
#                      ont           = "BP",
#                      qvalueCutoff  = 0.05,
#                      readable      = TRUE)
# head(ego_1200)
# 
# View(ego_1200@result)
# 
# barplot(ego_1200, showCategory = 20)


# ego_60_gbs<- enrichGO(gene          = rownames(res01_60_gbs),
#                      keyType = "ENSEMBL",
#                      universe   = rownames(cts),
#                      OrgDb         = org.Hs.eg.db,
#                      ont           = "BP",
#                      qvalueCutoff  = 0.05,
#                      readable      = TRUE)
# head(ego_1200)
# 
# View(ego_1200@result)
# 
# barplot(ego_60_gbs, showCategory = 20)
# 
# ego_1200_gbs <- enrichGO(gene          = rownames(res01_1200_gbs),
#                      keyType = "ENSEMBL",
#                      universe   = rownames(cts),
#                      OrgDb         = org.Hs.eg.db,
#                      ont           = "BP",
#                      qvalueCutoff  = 0.05,
#                      readable      = TRUE)
# head(ego_1200)
# 
# View(ego_1200@result)
# 
# barplot(ego_1200_gbs, showCategory = 20)


library(EnhancedVolcano)
p <- EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     y = 'pvalue',
                     xlim = c(-5, 5))

p + scale_y_continuous(
  name = bquote(~-Log[10]~italic(P)),
  limits=c(0, 50),
  breaks=c(0, 33, 40, 49))

library(GOplot)

# Generate the plotting object
genes <- as.data.frame(res01)
rownames(genes) <- NULL
#genes_names <- rownames(genes)
genes[genes == ""] = NA
sum(is.na(genes$hgnc_symbol))
genes <- na.omit(genes)
sum(is.na(genes$hgnc_symbol))
genes$ID <- genes$hgnc_symbol
genes <- genes[,c("log2FoldChange", "ID")]
colnames(genes) <- NULL
colnames(genes) <- c("logFC", "ID")
genes <- genes[,c("ID","logFC")]
terms <- as.data.frame(ego@result)
terms$geneID <- gsub("/", "," , terms$geneID)
terms <- terms[, c("ID", "Description", "p.adjust", "geneID")]
terms$category <- c(rep("BP", nrow(terms)))
colnames(terms) <- c("ID", "Term", "adj_pval","Genes", "category")
terms <- subset(terms, adj_pval < 0.05)
terms <- terms[, c("category","ID", "Term", "adj_pval","Genes")]
circ <- circle_dat(terms, genes)
GOCircle(circ)

process <- terms$Term
process <- process[c(1,2,6,8,7)]
chord <- chord_dat(circ,genes,process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)


# Generate the plotting object
genes_1200 <- as.data.frame(res01_1200)
rownames(genes_1200) <- NULL
#genes_names <- rownames(genes)
genes_1200[genes_1200 == ""] = NA
sum(is.na(genes_1200$hgnc_symbol))
genes_1200 <- na.omit(genes_1200)
sum(is.na(genes_1200$hgnc_symbol))
genes_1200$ID <- genes_1200$hgnc_symbol
genes_1200 <- genes_1200[,c("log2FoldChange", "ID")]
colnames(genes_1200) <- NULL
colnames(genes_1200) <- c("logFC", "ID")
genes_1200 <- genes_1200[,c("ID","logFC")]
terms_1200 <- as.data.frame(ego_1200@result)
terms_1200$geneID <- gsub("/", "," , terms_1200$geneID)
terms_1200 <- terms_1200[, c("ID", "Description", "p.adjust", "geneID")]
terms_1200$category <- c(rep("BP", nrow(terms)))
colnames(terms_1200) <- c("ID", "Term", "adj_pval","Genes", "category")
terms_1200 <- subset(terms_1200, adj_pval < 0.05)
terms_1200 <- terms_1200[, c("category","ID", "Term", "adj_pval","Genes")]
circ_1200 <- circle_dat(terms_1200, genes_1200)
GOCircle(circ_1200)

process_1200 <- terms_1200$Term
process_1200 <- process_1200[c(1,2,6,8,7)]
chord_1200 <- chord_dat(circ_1200,genes_1200,process_1200)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)





genes_60 <- as.data.frame(res01_60)
rownames(genes_60) <- NULL
#genes_names <- rownames(genes)
genes_60[genes_60 == ""] = NA
sum(is.na(genes_60$hgnc_symbol))
genes_60 <- na.omit(genes_60)
sum(is.na(genes_60$hgnc_symbol))
genes_60$ID <- genes_60$hgnc_symbol
genes_60 <- genes_60[,c("log2FoldChange", "ID")]
colnames(genes_60) <- NULL
colnames(genes_60) <- c("logFC", "ID")
genes_60 <- genes_60[,c("ID","logFC")]
terms_60 <- as.data.frame(ego_60@result)
terms_60$geneID <- gsub("/", "," , terms_60$geneID)
terms_60 <- terms_60[, c("ID", "Description", "p.adjust", "geneID")]
terms_60$category <- c(rep("BP", nrow(terms_60)))
colnames(terms_60) <- c("ID", "Term", "adj_pval","Genes", "category")
terms_60 <- subset(terms_60, adj_pval < 0.05)
terms_60 <- terms_60[, c("category","ID", "Term", "adj_pval","Genes")]
circ_60 <- circle_dat(terms_60, genes_60)
process_60 <- terms_60$Term
#process <- process[c(1,2,3,4,5,6,7)]
chord_60 <- chord_dat(circ_60,genes_60,process_60)
GOChord(chord_60, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)





library("pheatmap")
# this gives log2(n + 1)
#ntd <- normTransform(dds)
mat <- assay(rld)

topVarGenes <- order(rowVars(assay(rld)),decreasing=TRUE)[1:10000]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, annotation_col=df)
pheatmap(mat, cluster_rows=T, show_rownames = F, show_colnames = F,
         cluster_cols=T, annotation_col=df, annotation_legend = T)
idx <- rownames(res)[ which(res$padj < 0.001) ]
mat[ idx, ]
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=T)[1:1000]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(mat[select,], cluster_rows=T, show_rownames = F, show_colnames = F,
         cluster_cols=T, annotation_col=df, annotation_legend = T)
# df <- as.data.frame(colData(dds)[,c("condition","days_to_walk")])
# pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=F,
#          cluster_cols=FALSE, annotation_col=df, annotation_legend = T)
#########################################UPREGULATED############################################################################
# library(org.Hs.eg.db)
# library(clusterProfiler)
# egoup <- enrichGO(gene          = rownames(upres01),
#                 keyType = "ENSEMBL",
#                 universe   = rownames(res),
#                 OrgDb         = org.Hs.eg.db,
#                 ont           = "BP",
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)
# head(egoup)
# View(egoup@result)
#egores <- ego@result
#egores <- subset(egores, egores$p.adjust <0.05)
#emapplot(ego)
#dotplot(ego, showCategory = 20)
#barplot(egoup, showCategory = 20)
# DICTENSEMBL<-toTable(org.Hs.egENSEMBL)
# DICTENSEMBL1 <- subset(DICTENSEMBL$gene_id, DICTENSEMBL$ensembl_id %in% rownames(res01))
# kk <- enrichKEGG(gene         = DICTENSEMBL1,
#                  organism     = 'hsa',
#                  qvalueCutoff = 0.1
#                  )
# head(kk)
# View(kk@result)
# barplot(kk, drop=FALSE, showCategory=10)
# emapplot(kk)






# #res01<- subset(res, padj < 0.05)
# library(topGO)
# library(biomaRt)
# universe <- rownames(cts)
# genesOfInterest <- rownames(res01)
# species <- "Homo sapiens"
# ontology = "biological process"
# algorithm = "classic"
# statistic = "fisher"
# pValue = 0.05
# adjustMethod = "BH"
# 
# enriquece <- function(universe, genesOfInterest, species, ontology = "biological process",
#                       algorithm = "classic", statistic = "fisher", pValue = 0.05,
#                       adjustMethod = "BH"){
#   ontology <- tolower(ontology)
#   ontology <- gsub(" ", "_", ontology)
#   species <- tolower(species)
#   species <- gsub("^([[:alpha:]]).* ", "\\1", species)
#   species <- paste0(species, "_gene_ensembl")
#   ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = species)
#   GO <- biomaRt::getBM(filters = "ensembl_gene_id",
#                        attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
#                        values = universe, mart = ensembl)
#   GO[GO == ""] <- NA
#   GO <- stats::na.omit(GO)
#   GO <- GO[which(GO$namespace_1003 == ontology), c(1, 2)]
#   rm(ensembl)
#   gene2GO <- split(GO$go_id, GO$ensembl_gene_id)
#   gene2GO <- lapply(gene2GO, unique)
#   rm(species, GO)
#   ontology <- toupper(gsub("^([[:alpha:]]).*\\_([[:alpha:]]).*$", "\\1\\2", ontology))
#   geneList <- factor(as.integer(universe %in% genesOfInterest))
#   names(geneList) <- universe
#   myGOdata <- suppressMessages(methods::new("topGOdata",
#                                             ontology = ontology, allGenes = geneList,
#                                             annot = topGO::annFUN.gene2GO,
#                                             gene2GO = gene2GO))
#   result <- topGO::runTest(myGOdata, algorithm = algorithm,
#                            statistic = statistic)
#   result <- topGO::GenTable(myGOdata, result,
#                             topNodes = length(result@score))
#   colnames(result)[6] <- "pValue"
#   if(any(TRUE %in% grepl("^<", result[, "pValue"]))){
#     result$pValue <- gsub("^<", "", result$pValue)
#   }
#   result$pValue <- as.numeric(result$pValue)
#   result$pAdj <- stats::p.adjust(result[, "pValue"], method = adjustMethod)
#   result <- result[result$pAdj <= pValue, ]
#   rownames(result) <- NULL
#   return(result)
# }
# enrich <- enriquece(rownames(cts), rownames(res01), "Homo sapiens", "biological process")
# idx05 <- enrich$pAdj < 0.05
# enrich05 <- enrich[idx05,]
##################################################transcriptogrammer###############################################################
# library(data.table)
# c1 <- fread("counts_npadcoding.txt", stringsAsFactors = F)
# c2 <- fread("counts_iowa_fastp.txt", stringsAsFactors = F)
# counts <- merge(c1, c2, by = "Geneid")
# counts <- as.data.frame(counts)
# rownames(counts) <- counts$Geneid
# counts <- counts[,-c(1:6,37:40,  97:101)]
# #c1 <- cbind(c1, c2[,-c(1:6)])
# #rm(c2)
# #c1 <- c1[,-c(2:6)]
# #counts <- as.data.frame(c1)
# #rm(c1)
# library(limma)
# library(edgeR)
# 
# #counts <- counts[, c(-1,-37, -38, -39, -40)]
# cpms <- cpm(counts)
# #keep <- filterByExpr(dge, design)
# keep <- rowSums(cpms >= 1) >= 1
# counts <- counts[keep,]
# nrow(counts)
# groups <- c(rep("REC", 4),#10##Demyelinating
#             rep("GBS", 4),#11##Demyelinating
#             rep("REC", 4),#12##Miller_Fisher
#             rep("GBS", 4),#13##Axonal
#             rep("REC", 4),#14##Axonal
#             rep("GBS", 4),#15##Demyelinating
#             rep("GBS", 4),#16##Miller_Fisher
#             rep("REC", 4),#17##Miller_Fisher
#             rep("GBS", 4),#18##Axonal
#             rep("GBS", 4),#1##Demyelinating
#             rep("REC", 4),#20##Demyelinating
#             rep("GBS", 4),#21##Demyelinating
#             rep("REC", 4),#22##Demyelinating
#             rep("REC", 4),#23##Demyelinating
#             rep("GBS", 4),#24##Demyelinating
#             rep("REC", 4),#2##Demyelinating
#             rep("GBS", 4),#3##Demyelinating
#             rep("REC", 4),#4##Demyelinating
#             rep("GBS", 4),#5##Axonal
#             rep("REC", 4),#6##Axonal
#             rep("GBS", 4),#7##Axonal
#             rep("REC", 4),#8##Axonal
#             rep("GBS", 4),#9##Miller_Fisher
#             rep("CTR", 3),#25##Control
#             rep("CTR", 3),#26##Control
#             rep("CTR", 3),#27##Control
#             rep("CTR", 3),#28##Control
#             rep("CTR", 3),#29##Control
#             rep("REC", 3),#30##Unknown
#             rep("REC", 3),#31##Unknown
#             rep("REC", 3),#32##Unknown
#             rep("REC", 3),#33##Unknown
#             rep("REC", 3),#34##Unknown
#             rep("REC", 3),#35##Unknown
#             rep("REC", 3),#36##Unknown
#             rep("REC", 3),#37##Unknown
#             rep("REC", 3),#38##Axonal
#             rep("GBS", 3),#39##Axonal
#             rep("GBS", 3),#40##Demyelinating
#             rep("REC", 3),#41##Demyelinating
#             rep("GBS", 3),#42##Axonal
#             rep("GBS", 3),#43##Axonal
#             rep("REC", 3),#44##Axonal
#             rep("REC", 3),#45##Demyelinating
#             rep("GBS", 3),#46##Demyelinating
#             rep("REC", 3),#47##Axonal
#             rep("GBS", 3)#48##Miller_Fisher
# )
# dge <- DGEList(counts, group = groups)
# dge <- calcNormFactors(dge)
# # plotMDS(dge, labels = colnames(counts), col = c("red", "darkgreen",
# #                                                 "blue")[factor(groups)])
# groups <- as.factor(groups)
# design <- model.matrix(~0 + groups)
# v <- voom(dge, design, plot=FALSE)
# logCPM <- v$E
# library(biomaRt)
# ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# dict <- getBM(filters = "ensembl_gene_id",
#               attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
#               values = rownames(logCPM), mart = ensembl)
# dict[dict==""] <- NA
# dict <- na.omit(dict)
# dict$ensembl_peptide_id <- paste0("9606.", dict$ensembl_peptide_id)
# links <- fread("9606.protein.links.v11.0.txt", stringsAsFactors = F)
# links <- links[links$combined_score >= 700, c(1, 2)]
# rm(counts, cpms, design, dge, ensembl, v, keep)
# library(transcriptogramer)
# library(RTN)
# library(igraph)
# library(RedeR)
# logCPM <- as.data.frame(logCPM)
# logCPM <- logCPM[, groups != "REC"]
# t <- transcriptogramPreprocess(association = links, ordering = Hs700)
# t <- transcriptogramStep1(object = t, expression = logCPM,
#                           dictionary = dict, nCores = T)
# radius(object = t) <- 80
# t <- transcriptogramStep2(object = t, nCores = T)
# groups <- as.character(groups)
# groups <- groups[groups!="REC"]
# groups[groups=="CTR"] <- T
# groups[groups=="GBS"] <- F
# groups <- as.logical(groups)
# t <- differentiallyExpressed(object = t, levels = groups, pValue = 0.001, foldChange = 0.13, species = "Homo sapiens")
# DE <- DE(t)
# terms <- clusterEnrichment(t, species = "Homo sapiens", pValue = 0.0001, nCores = T)
# termos <- terms@Terms
# clustergenes <- t@DE
# rdp <- clusterVisualization(object = t)

# ##############################################AMOSTRAS menos o outlier  sample 19  #############################################
# cts2 <- cts[,-c(37:40)]
# #library("car")
# #which.names(c("19.L5", "19.L6", "19.L7", "19.L8"), coldata)
# coldata2 <-coldata[-c(37:40), ]
# all(rownames(coldata2) %in% colnames(cts2))
# all(rownames(coldata2) == colnames(cts2))
# library(DESeq2)
# dds2 <- DESeqDataSetFromMatrix(countData = cts2,
#                                colData = coldata2,
#                                design = ~ condition)
# 
# featureData2 <- data.frame(gene=rownames(cts2))
# mcols(dds2) <- DataFrame(mcols(dds2), featureData2)
# mcols(dds2)
# keep2 <- rowSums(counts(dds2)) >= 10
# dds2 <- dds2[keep2,]
# dds2$condition <- factor(dds2$condition, levels = c("REC","GBS", "CTR"))
# dds2$condition <- relevel(dds2$condition, ref = "CTR")
# dds2$condition <- droplevels(dds2$condition)
# dds2 <- DESeq(dds2)
# res2 <- results(dds2)
# res2
# res2 <- results(dds2, name="condition_GBS_vs_CTR")
# res2 <- results(dds2, contrast=c("condition","GBS","CTR"))
# nrow(res2)
# head( res2[ order( res2$log2FoldChange ), ] )
# plotMA( res2, ylim = c(-4, 4) )
# res02 <- subset(res2, padj < 0.001 & abs(log2FoldChange)> 0.6)
# nrow(res02)
# #res01 <- subset(res, padj < 0.05)
# plotMA( res02, ylim = c(-4, 4) )
# nrow(res02)
# head( res02[ order( res01$log2FoldChange ), ] )
# #plotDispEsts( dds, ylim = c(1e-6, 1e3) )
# res02 <- as.data.frame(res02)
# upres02 <- subset(res02, res02$log2FoldChange > 0)
# downres02 <- subset(res02, res02$log2FoldChange < 0)
# res2$ensembl <- sapply( strsplit( rownames(res2), split="\\+" ), "[", 1 )
# library( "biomaRt" )
# ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
# genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
#                   filters = "ensembl_gene_id",
#                   values = res$ensembl,
#                   mart = ensembl )
# idx2 <- match( res2$ensembl, genemap$ensembl_gene_id )
# res2$entrez <- genemap$entrezgene[ idx2 ]
# res2$hgnc_symbol <- genemap$hgnc_symbol[ idx2 ]
# #res02 <- subset(res2, padj < 0.05 & abs(log2FoldChange)>1)
# nrow(res02)
# library("ggplot2")
# library("RColorBrewer")
# vsd2<- vst(dds2, blind=FALSE)
# 
# pcaData2 <- plotPCA(vsd2, intgroup=c("condition", "type"), ntop= 21000, returnData=TRUE)
# percentVar2 <- round(100 * attr(pcaData2, "percentVar"))
# ggplot(pcaData2, aes(PC1, PC2, color=condition, shape=type)) +
#   geom_point(size=7) +
#   xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar2[2],"% variance")) + geom_text(aes_string(x = "PC1", y = "PC2", label = "name"), color = "black", nudge_x = 2, nudge_y = -2, size = 3) +
#   coord_fixed()
# 
# library(org.Hs.eg.db)
# library(clusterProfiler)
# ego2 <- enrichGO(gene          = rownames(res02),
#                  keyType = "ENSEMBL",
#                  universe   = rownames(cts2),
#                  OrgDb         = org.Hs.eg.db,
#                  ont           = "BP",
#                  qvalueCutoff  = 0.05,
#                  readable      = TRUE)
# head(ego2)
# View(ego2@result)
# #egores2 <- subset(egores, egores2$p.adjust <0.05)
# #emapplot(ego2)
# dotplot(ego2, showCategory = 20)
# barplot(ego2, showCategory =40)
# 
# DICTENSEMBL<-toTable(org.Hs.egENSEMBL)
# DICTENSEMBL2 <- subset(DICTENSEMBL$gene_id, DICTENSEMBL$ensembl_id %in% rownames(res02))
# kk2 <- enrichKEGG(gene         = DICTENSEMBL2,
#                   organism     = 'hsa',
#                   pvalueCutoff = 0.05)
# head(kk2)
# View(kk2@result)
# barplot(kk2, drop=FALSE, showCategory=40)
# emapplot(kk)
# #######################################################TOPGO_res02####################################
# library(topGO)
# library(biomaRt)
# universe <- rownames(cts2)
# genesOfInterest <- rownames(res02)
# species <- "Homo sapiens"
# ontology = "biological process"
# algorithm = "classic"
# statistic = "fisher"
# pValue = 0.05
# adjustMethod = "BH"
# 
# enriquece <- function(universe, genesOfInterest, species, ontology = "biological process",
#                       algorithm = "classic", statistic = "fisher", pValue = 0.05,
#                       adjustMethod = "BH"){
#   ontology <- tolower(ontology)
#   ontology <- gsub(" ", "_", ontology)
#   species <- tolower(species)
#   species <- gsub("^([[:alpha:]]).* ", "\\1", species)
#   species <- paste0(species, "_gene_ensembl")
#   ensembl <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL", dataset = species)
#   GO <- biomaRt::getBM(filters = "ensembl_gene_id",
#                        attributes = c("ensembl_gene_id", "go_id", "namespace_1003"),
#                        values = universe, mart = ensembl)
#   GO[GO == ""] <- NA
#   GO <- stats::na.omit(GO)
#   GO <- GO[which(GO$namespace_1003 == ontology), c(1, 2)]
#   rm(ensembl)
#   gene2GO <- split(GO$go_id, GO$ensembl_gene_id)
#   gene2GO <- lapply(gene2GO, unique)
#   rm(species, GO)
#   ontology <- toupper(gsub("^([[:alpha:]]).*\\_([[:alpha:]]).*$", "\\1\\2", ontology))
#   geneList <- factor(as.integer(universe %in% genesOfInterest))
#   names(geneList) <- universe
#   myGOdata <- suppressMessages(methods::new("topGOdata",
#                                             ontology = ontology, allGenes = geneList,
#                                             annot = topGO::annFUN.gene2GO,
#                                             gene2GO = gene2GO))
#   result <- topGO::runTest(myGOdata, algorithm = algorithm,
#                            statistic = statistic)
#   result <- topGO::GenTable(myGOdata, result,
#                             topNodes = length(result@score))
#   colnames(result)[6] <- "pValue"
#   if(any(TRUE %in% grepl("^<", result[, "pValue"]))){
#     result$pValue <- gsub("^<", "", result$pValue)
#   }
#   result$pValue <- as.numeric(result$pValue)
#   result$pAdj <- stats::p.adjust(result[, "pValue"], method = adjustMethod)
#   result <- result[result$pAdj <= pValue, ]
#   rownames(result) <- NULL
#   return(result)
# }
# enrich2 <- enriquece(rownames(cts2), rownames(res02), "Homo sapiens", "biological process")
# 
# ##################################################transcriptogrammer_res2###############################################################
# library(data.table)
# c1s <- fread("counts_npadcoding.txt", stringsAsFactors = F)
# c2s <- fread("counts_coding.txt", stringsAsFactors = F)
# c1s <- cbind(c1s, c2s[,-c(1:6)])
# rm(c2s)
# c1s <- c1s[,-c(2:6)]
# c1s <- c1s[,-c(38:41)]
# counts <- as.data.frame(c1s)
# rm(c1s)
# library(limma)
# library(edgeR)
# rownames(counts) <- counts$Geneid
# counts <- counts[, -1]
# cpms <- cpm(counts)
# #keep <- filterByExpr(dge, design)
# keep <- rowSums(cpms >= 1) >= 1
# counts <- counts[keep,]
# nrow(counts)
# groups <- c(rep("REC", 4),#10##Demyelinating
#             rep("GBS", 4),#11##Demyelinating
#             rep("REC", 4),#12##Miller_Fisher
#             rep("GBS", 4),#13##Axonal
#             rep("REC", 4),#14##Axonal
#             rep("GBS", 4),#15##Demyelinating
#             rep("GBS", 4),#16##Miller_Fisher
#             rep("REC", 4),#17##Miller_Fisher
#             rep("GBS", 4),#18##Axonal
#             rep("GBS", 4),#1##Demyelinating
#             rep("REC", 4),#20##Demyelinating
#             rep("GBS", 4),#21##Demyelinating
#             rep("REC", 4),#22##Demyelinating
#             rep("REC", 4),#23##Demyelinating
#             rep("GBS", 4),#24##Demyelinating
#             rep("REC", 4),#2##Demyelinating
#             rep("GBS", 4),#3##Demyelinating
#             rep("REC", 4),#4##Demyelinating
#             rep("GBS", 4),#5##Axonal
#             rep("REC", 4),#6##Axonal
#             rep("GBS", 4),#7##Axonal
#             rep("REC", 4),#8##Axonal
#             rep("GBS", 4),#9##Miller_Fisher
#             rep("CTR", 3),#25##Control
#             rep("CTR", 3),#26##Control
#             rep("CTR", 3),#27##Control
#             rep("CTR", 3),#28##Control
#             rep("CTR", 3),#29##Control
#             rep("REC", 3),#30##Unknown
#             rep("REC", 3),#31##Unknown
#             rep("REC", 3),#32##Unknown
#             rep("REC", 3),#33##Unknown
#             rep("REC", 3),#34##Unknown
#             rep("REC", 3),#35##Unknown
#             rep("REC", 3),#36##Unknown
#             rep("REC", 3),#37##Unknown
#             rep("REC", 3),#38##Axonal
#             rep("GBS", 3),#39##Axonal
#             rep("GBS", 3),#40##Demyelinating
#             rep("REC", 3),#41##Demyelinating
#             rep("GBS", 3),#42##Axonal
#             rep("GBS", 3),#43##Axonal
#             rep("REC", 3),#44##Axonal
#             rep("REC", 3),#45##Demyelinating
#             rep("GBS", 3),#46##Demyelinating
#             rep("REC", 3),#47##Axonal
#             rep("GBS", 3)#48##Miller_Fisher
# )
# dge <- DGEList(counts, group = groups)
# dge <- calcNormFactors(dge)
# # plotMDS(dge, labels = colnames(counts), col = c("red", "darkgreen",
# #                                                 "blue")[factor(groups)])
# groups <- as.factor(groups)
# design <- model.matrix(~0 + groups)
# v <- voom(dge, design, plot=FALSE)
# logCPM <- v$E
# library(biomaRt)
# ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# dict <- getBM(filters = "ensembl_gene_id",
#               attributes = c("ensembl_peptide_id", "ensembl_gene_id"),
#               values = rownames(logCPM), mart = ensembl)
# dict[dict==""] <- NA
# dict <- na.omit(dict)
# dict$ensembl_peptide_id <- paste0("9606.", dict$ensembl_peptide_id)
# links <- fread("9606.protein.links.v11.0.txt", stringsAsFactors = F)
# links <- links[links$combined_score >= 700, c(1, 2)]
# rm(counts, cpms, design, dge, ensembl, v, keep)
# library(transcriptogramer)
# library(RTN)
# library(igraph)