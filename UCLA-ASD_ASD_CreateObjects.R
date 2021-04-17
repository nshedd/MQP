library(dplyr)
library(Seurat)
library(ggplot2)

print('1BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/1BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW1 = CreateSeuratObject(counts = expression_matrix, project='1BW')
BW1[["percent.mt"]] <- PercentageFeatureSet(BW1, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW1$nCount_RNA)>=minCov){
  countLOW=min(BW1$nCount_RNA)
}else{
  countLOW=quantile(BW1$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW1$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW1$nFeature_RNA, prob=0.01)
BW1 <- subset(BW1, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW1$Sample <- "BW1"


print('3BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/3BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW3 = CreateSeuratObject(counts = expression_matrix, project='3BW')
BW3[["percent.mt"]] <- PercentageFeatureSet(BW3, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW3$nCount_RNA)>=minCov){
  countLOW=min(BW3$nCount_RNA)
}else{
  countLOW=quantile(BW3$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW3$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW3$nFeature_RNA, prob=0.01)
BW3 <- subset(BW3, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW3$Sample <- "BW3"


print('4BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/4BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW4 = CreateSeuratObject(counts = expression_matrix, project='4BW')
BW4[["percent.mt"]] <- PercentageFeatureSet(BW4, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW4$nCount_RNA)>=minCov){
  countLOW=min(BW4$nCount_RNA)
}else{
  countLOW=quantile(BW4$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW4$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW4$nFeature_RNA, prob=0.01)
BW4 <- subset(BW4, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW4$Sample <- "BW4"


print('7BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/7BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW7 = CreateSeuratObject(counts = expression_matrix, project='7BW')
BW7[["percent.mt"]] <- PercentageFeatureSet(BW7, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW7$nCount_RNA)>=minCov){
  countLOW=min(BW7$nCount_RNA)
}else{
  countLOW=quantile(BW7$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW7$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW7$nFeature_RNA, prob=0.01)
BW7 <- subset(BW7, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW7$Sample <- "BW7"

print('9BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/9BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW9 = CreateSeuratObject(counts = expression_matrix, project='9BW')
BW9[["percent.mt"]] <- PercentageFeatureSet(BW9, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW9$nCount_RNA)>=minCov){
  countLOW=min(BW9$nCount_RNA)
}else{
  countLOW=quantile(BW9$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW9$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW9$nFeature_RNA, prob=0.01)
BW9 <- subset(BW9, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW9$Sample <- "BW9"


print('10BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/10BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW10 = CreateSeuratObject(counts = expression_matrix, project='10BW')
BW10[["percent.mt"]] <- PercentageFeatureSet(BW10, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW10$nCount_RNA)>=minCov){
  countLOW=min(BW10$nCount_RNA)
}else{
  countLOW=quantile(BW10$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW10$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW10$nFeature_RNA, prob=0.01)
BW10 <- subset(BW10, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW10$Sample <- "BW10"


print('14BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/14BW_S17/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW14 = CreateSeuratObject(counts = expression_matrix, project='14BW')
BW14[["percent.mt"]] <- PercentageFeatureSet(BW14, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW14$nCount_RNA)>=minCov){
  countLOW=min(BW14$nCount_RNA)
}else{
  countLOW=quantile(BW14$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW14$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW14$nFeature_RNA, prob=0.01)
BW14 <- subset(BW14, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW14$Sample <- "BW14"


print('16BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/16BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW16 = CreateSeuratObject(counts = expression_matrix, project='16BW')
BW16[["percent.mt"]] <- PercentageFeatureSet(BW16, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW16$nCount_RNA)>=minCov){
  countLOW=min(BW16$nCount_RNA)
}else{
  countLOW=quantile(BW16$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW16$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW16$nFeature_RNA, prob=0.01)
BW16 <- subset(BW16, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW16$Sample <- "BW16"


print('17BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/17BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW17 = CreateSeuratObject(counts = expression_matrix, project='17BW')
BW17[["percent.mt"]] <- PercentageFeatureSet(BW17, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW17$nCount_RNA)>=minCov){
  countLOW=min(BW17$nCount_RNA)
}else{
  countLOW=quantile(BW17$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW17$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW17$nFeature_RNA, prob=0.01)
BW17 <- subset(BW17, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW17$Sample <- "BW17"

print('18BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/18BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW18 = CreateSeuratObject(counts = expression_matrix, project='18BW')
BW18[["percent.mt"]] <- PercentageFeatureSet(BW18, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW18$nCount_RNA)>=minCov){
  countLOW=min(BW18$nCount_RNA)
}else{
  countLOW=quantile(BW18$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW18$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW18$nFeature_RNA, prob=0.01)
BW18 <- subset(BW18, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW18$Sample <- "BW18"

print('19BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/19BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW19 = CreateSeuratObject(counts = expression_matrix, project='19BW')
BW19[["percent.mt"]] <- PercentageFeatureSet(BW19, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW19$nCount_RNA)>=minCov){
  countLOW=min(BW19$nCount_RNA)
}else{
  countLOW=quantile(BW19$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW19$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW19$nFeature_RNA, prob=0.01)
BW19 <- subset(BW19, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW19$Sample <- "BW19"

print('21BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/21BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW21 = CreateSeuratObject(counts = expression_matrix, project='21BW')
BW21[["percent.mt"]] <- PercentageFeatureSet(BW21, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW21$nCount_RNA)>=minCov){
  countLOW=min(BW21$nCount_RNA)
}else{
  countLOW=quantile(BW21$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW21$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW21$nFeature_RNA, prob=0.01)
BW21 <- subset(BW21, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW21$Sample <- "BW21"

print('22BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/22BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW22 = CreateSeuratObject(counts = expression_matrix, project='22BW')
BW22[["percent.mt"]] <- PercentageFeatureSet(BW22, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW22$nCount_RNA)>=minCov){
  countLOW=min(BW22$nCount_RNA)
}else{
  countLOW=quantile(BW22$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW22$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW22$nFeature_RNA, prob=0.01)
BW22 <- subset(BW22, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW22$Sample <- "BW22"

print('24BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/24BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW24 = CreateSeuratObject(counts = expression_matrix, project='24BW')
BW24[["percent.mt"]] <- PercentageFeatureSet(BW24, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW24$nCount_RNA)>=minCov){
  countLOW=min(BW24$nCount_RNA)
}else{
  countLOW=quantile(BW24$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW24$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW24$nFeature_RNA, prob=0.01)
BW24 <- subset(BW24, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW24$Sample <- "BW24"


print('26BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/26BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW26 = CreateSeuratObject(counts = expression_matrix, project='26BW')
BW26[["percent.mt"]] <- PercentageFeatureSet(BW26, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW26$nCount_RNA)>=minCov){
  countLOW=min(BW26$nCount_RNA)
}else{
  countLOW=quantile(BW26$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW26$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW26$nFeature_RNA, prob=0.01)
BW26 <- subset(BW26, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW26$Sample <- "BW26"

print('27BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/27BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW27 = CreateSeuratObject(counts = expression_matrix, project='27BW')
BW27[["percent.mt"]] <- PercentageFeatureSet(BW27, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW27$nCount_RNA)>=minCov){
  countLOW=min(BW27$nCount_RNA)
}else{
  countLOW=quantile(BW27$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW27$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW27$nFeature_RNA, prob=0.01)
BW27 <- subset(BW27, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW27$Sample <- "BW27"

print('28BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/28BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW28 = CreateSeuratObject(counts = expression_matrix, project='28BW')
BW28[["percent.mt"]] <- PercentageFeatureSet(BW28, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW28$nCount_RNA)>=minCov){
  countLOW=min(BW28$nCount_RNA)
}else{
  countLOW=quantile(BW28$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW28$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW28$nFeature_RNA, prob=0.01)
BW28 <- subset(BW28, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW28$Sample <- "BW28"

print('29BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/29BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW29 = CreateSeuratObject(counts = expression_matrix, project='29BW')
BW29[["percent.mt"]] <- PercentageFeatureSet(BW29, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW29$nCount_RNA)>=minCov){
  countLOW=min(BW29$nCount_RNA)
}else{
  countLOW=quantile(BW29$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW29$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW29$nFeature_RNA, prob=0.01)
BW29 <- subset(BW29, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW29$Sample <- "BW29"

print('31BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/31BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW31 = CreateSeuratObject(counts = expression_matrix, project='31BW')
BW31[["percent.mt"]] <- PercentageFeatureSet(BW31, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW31$nCount_RNA)>=minCov){
  countLOW=min(BW31$nCount_RNA)
}else{
  countLOW=quantile(BW31$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW31$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW31$nFeature_RNA, prob=0.01)
BW31 <- subset(BW31, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW31$Sample <- "BW31"

print('32BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/32BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW32 = CreateSeuratObject(counts = expression_matrix, project='32BW')
BW32[["percent.mt"]] <- PercentageFeatureSet(BW32, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW32$nCount_RNA)>=minCov){
  countLOW=min(BW32$nCount_RNA)
}else{
  countLOW=quantile(BW32$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW32$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW32$nFeature_RNA, prob=0.01)
BW32 <- subset(BW32, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW32$Sample <- "BW32"

print('34BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/34BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW34 = CreateSeuratObject(counts = expression_matrix, project='34BW')
BW34[["percent.mt"]] <- PercentageFeatureSet(BW34, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW34$nCount_RNA)>=minCov){
  countLOW=min(BW34$nCount_RNA)
}else{
  countLOW=quantile(BW34$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW34$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW34$nFeature_RNA, prob=0.01)
BW34 <- subset(BW34, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW34$Sample <- "BW34"

print('37BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/37BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW37 = CreateSeuratObject(counts = expression_matrix, project='37BW')
BW37[["percent.mt"]] <- PercentageFeatureSet(BW37, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW37$nCount_RNA)>=minCov){
  countLOW=min(BW37$nCount_RNA)
}else{
  countLOW=quantile(BW37$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW37$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW37$nFeature_RNA, prob=0.01)
BW37 <- subset(BW37, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW37$Sample <- "BW37"

print('38BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/38BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW38 = CreateSeuratObject(counts = expression_matrix, project='38BW')
BW38[["percent.mt"]] <- PercentageFeatureSet(BW38, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW38$nCount_RNA)>=minCov){
  countLOW=min(BW38$nCount_RNA)
}else{
  countLOW=quantile(BW38$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW38$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW38$nFeature_RNA, prob=0.01)
BW38 <- subset(BW38, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW38$Sample <- "BW38"

print('39BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/39BW_S9/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW39 = CreateSeuratObject(counts = expression_matrix, project='39BW')
BW39[["percent.mt"]] <- PercentageFeatureSet(BW39, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW39$nCount_RNA)>=minCov){
  countLOW=min(BW39$nCount_RNA)
}else{
  countLOW=quantile(BW39$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW39$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW39$nFeature_RNA, prob=0.01)
BW39 <- subset(BW39, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW39$Sample <- "BW39"

print('42BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/42BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW42 = CreateSeuratObject(counts = expression_matrix, project='42BW')
BW42[["percent.mt"]] <- PercentageFeatureSet(BW42, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW42$nCount_RNA)>=minCov){
  countLOW=min(BW42$nCount_RNA)
}else{
  countLOW=quantile(BW42$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW42$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW42$nFeature_RNA, prob=0.01)
BW42 <- subset(BW42, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW42$Sample <- "BW42"

print('44BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/44BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW44 = CreateSeuratObject(counts = expression_matrix, project='44BW')
BW44[["percent.mt"]] <- PercentageFeatureSet(BW44, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW44$nCount_RNA)>=minCov){
  countLOW=min(BW44$nCount_RNA)
}else{
  countLOW=quantile(BW44$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW44$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW44$nFeature_RNA, prob=0.01)
BW44 <- subset(BW44, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW44$Sample <- "BW44"

print('45BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/45BW_S15/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW45 = CreateSeuratObject(counts = expression_matrix, project='45BW')
BW45[["percent.mt"]] <- PercentageFeatureSet(BW45, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW45$nCount_RNA)>=minCov){
  countLOW=min(BW45$nCount_RNA)
}else{
  countLOW=quantile(BW45$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW45$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW45$nFeature_RNA, prob=0.01)
BW45 <- subset(BW45, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW45$Sample <- "BW45"

print('48BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/48BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW48 = CreateSeuratObject(counts = expression_matrix, project='48BW')
BW48[["percent.mt"]] <- PercentageFeatureSet(BW48, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW48$nCount_RNA)>=minCov){
  countLOW=min(BW48$nCount_RNA)
}else{
  countLOW=quantile(BW48$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW48$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW48$nFeature_RNA, prob=0.01)
BW48 <- subset(BW48, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW48$Sample <- "BW48"

print('49BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/49BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW49 = CreateSeuratObject(counts = expression_matrix, project='49BW')
BW49[["percent.mt"]] <- PercentageFeatureSet(BW49, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW49$nCount_RNA)>=minCov){
  countLOW=min(BW49$nCount_RNA)
}else{
  countLOW=quantile(BW49$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW49$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW49$nFeature_RNA, prob=0.01)
BW49 <- subset(BW49, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW49$Sample <- "BW49"

print('50BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/50BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW50 = CreateSeuratObject(counts = expression_matrix, project='50BW')
BW50[["percent.mt"]] <- PercentageFeatureSet(BW50, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW50$nCount_RNA)>=minCov){
  countLOW=min(BW50$nCount_RNA)
}else{
  countLOW=quantile(BW50$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW50$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW50$nFeature_RNA, prob=0.01)
BW50 <- subset(BW50, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW50$Sample <- "BW50"

print('54BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/54BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW54 = CreateSeuratObject(counts = expression_matrix, project='54BW')
BW54[["percent.mt"]] <- PercentageFeatureSet(BW54, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW54$nCount_RNA)>=minCov){
  countLOW=min(BW54$nCount_RNA)
}else{
  countLOW=quantile(BW54$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW54$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW54$nFeature_RNA, prob=0.01)
BW54 <- subset(BW54, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW54$Sample <- "BW54"

print('59BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/59BW_S29/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW59 = CreateSeuratObject(counts = expression_matrix, project='59BW')
BW59[["percent.mt"]] <- PercentageFeatureSet(BW59, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW59$nCount_RNA)>=minCov){
  countLOW=min(BW59$nCount_RNA)
}else{
  countLOW=quantile(BW59$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW59$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW59$nFeature_RNA, prob=0.01)
BW59 <- subset(BW59, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW59$Sample <- "BW59"

print('60BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/60BW_S30/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW60 = CreateSeuratObject(counts = expression_matrix, project='60BW')
BW60[["percent.mt"]] <- PercentageFeatureSet(BW60, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW60$nCount_RNA)>=minCov){
  countLOW=min(BW60$nCount_RNA)
}else{
  countLOW=quantile(BW60$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW60$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW60$nFeature_RNA, prob=0.01)
BW60 <- subset(BW60, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW60$Sample <- "BW60"

print('65BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/65BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW65 = CreateSeuratObject(counts = expression_matrix, project='65BW')
BW65[["percent.mt"]] <- PercentageFeatureSet(BW65, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW65$nCount_RNA)>=minCov){
  countLOW=min(BW65$nCount_RNA)
}else{
  countLOW=quantile(BW65$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW65$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW65$nFeature_RNA, prob=0.01)
BW65 <- subset(BW65, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW65$Sample <- "BW65"

print('66BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/66BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW66 = CreateSeuratObject(counts = expression_matrix, project='66BW')
BW66[["percent.mt"]] <- PercentageFeatureSet(BW66, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW66$nCount_RNA)>=minCov){
  countLOW=min(BW66$nCount_RNA)
}else{
  countLOW=quantile(BW66$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW66$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW66$nFeature_RNA, prob=0.01)
BW66 <- subset(BW66, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW66$Sample <- "BW66"

print('67BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/67BW_S9/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW67 = CreateSeuratObject(counts = expression_matrix, project='67BW')
BW67[["percent.mt"]] <- PercentageFeatureSet(BW67, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW67$nCount_RNA)>=minCov){
  countLOW=min(BW67$nCount_RNA)
}else{
  countLOW=quantile(BW67$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW67$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW67$nFeature_RNA, prob=0.01)
BW67 <- subset(BW67, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW67$Sample <- "BW67"

print('68BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/68BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW68 = CreateSeuratObject(counts = expression_matrix, project='68BW')
BW68[["percent.mt"]] <- PercentageFeatureSet(BW68, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW68$nCount_RNA)>=minCov){
  countLOW=min(BW68$nCount_RNA)
}else{
  countLOW=quantile(BW68$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW68$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW68$nFeature_RNA, prob=0.01)
BW68 <- subset(BW68, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW68$Sample <- "BW68"

print('70BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/70BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW70 = CreateSeuratObject(counts = expression_matrix, project='70BW')
BW70[["percent.mt"]] <- PercentageFeatureSet(BW70, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW70$nCount_RNA)>=minCov){
  countLOW=min(BW70$nCount_RNA)
}else{
  countLOW=quantile(BW70$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW70$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW70$nFeature_RNA, prob=0.01)
BW70 <- subset(BW70, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW70$Sample <- "BW70"

print('72BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/72BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW72 = CreateSeuratObject(counts = expression_matrix, project='72BW')
BW72[["percent.mt"]] <- PercentageFeatureSet(BW72, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW72$nCount_RNA)>=minCov){
  countLOW=min(BW72$nCount_RNA)
}else{
  countLOW=quantile(BW72$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW72$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW72$nFeature_RNA, prob=0.01)
BW72 <- subset(BW72, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW72$Sample <- "BW72"

print('73BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/73BW_S15/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW73 = CreateSeuratObject(counts = expression_matrix, project='73BW')
BW73[["percent.mt"]] <- PercentageFeatureSet(BW73, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW73$nCount_RNA)>=minCov){
  countLOW=min(BW73$nCount_RNA)
}else{
  countLOW=quantile(BW73$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW73$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW73$nFeature_RNA, prob=0.01)
BW73 <- subset(BW73, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW73$Sample <- "BW73"

print('75BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/75BW_S17/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW75 = CreateSeuratObject(counts = expression_matrix, project='75BW')
BW75[["percent.mt"]] <- PercentageFeatureSet(BW75, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW75$nCount_RNA)>=minCov){
  countLOW=min(BW75$nCount_RNA)
}else{
  countLOW=quantile(BW75$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW75$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW75$nFeature_RNA, prob=0.01)
BW75 <- subset(BW75, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW75$Sample <- "BW75"

print('76BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/76BW_S18/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW76 = CreateSeuratObject(counts = expression_matrix, project='76BW')
BW76[["percent.mt"]] <- PercentageFeatureSet(BW76, pattern = "^MT-")
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW76$nCount_RNA)>=minCov){
  countLOW=min(BW76$nCount_RNA)
}else{
  countLOW=quantile(BW76$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW76$nCount_RNA, prob=0.95) 
featureLOW=quantile(BW76$nFeature_RNA, prob=0.01)
BW76 <- subset(BW76, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH & percent.mt < 5)
BW76$Sample <- "BW76"



CTL_BA4.6 <- merge(BW60, y=c(BW31,BW26),
                   add.cell.ids=c('BW60','BW31','BW26'),
                   project = "UCLA-ASD")
                   #3.7.5,  3.6.7,  3.6.5   3.5.4
                   #8751,   7152,   7076,   4829
                   #2593,   373,    1519,   18072
saveRDS(CTL_BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/data/ASD_BA4.6')

CTL_BA9 <- merge(BW66, y=c(BW34, BW54),
                 add.cell.ids=c('BW66','BW34','BW54'),
                 project = "UCLA-ASD")
                 #3.9.7,  3.6.5,  3.11.10
                 #10242,   8453,     6818

saveRDS(CTL_BA9, '/data/rusers/sheddn/UCLA-ASD/subset/data/ASD_BA9')

q()

ASD_BA4.6 <- merge(BW1, y=c(BW10,BW14,BW16,BW17,BW22,BW24,BW26,BW27,BW3,BW31,BW32,BW4,BW59,BW60,BW7,BW9),
                   add.cell.ids=c('BW1','BW10','BW14','BW16','BW17','BW22','BW24','BW26','BW27','BW3','BW31','BW32','BW4','BW59','BW60','BW7','BW9'),
                   project = "UCLA-ASD")

saveRDS(ASD_BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/ASD_BA4.6')

ASD_BA9 <- merge(BW18, y=c(BW19,BW21,BW28,BW29,BW34,BW37,BW38,BW39,BW42,BW44,BW45,BW48,BW49,BW50,BW54,BW65,BW66,BW67,BW68,BW70,BW72,BW73,BW75,BW76),
                 add.cell.ids=c('BW18','BW19','BW21','BW28','BW29','BW34','BW37','BW38','BW39','BW42','BW44','BW45','BW48','BW49','BW50','BW54','BW65','BW66','BW67','BW68','BW70','BW72','BW73','BW75','BW76'),
                 project = "UCLA-ASD")

saveRDS(ASD_BA9, '/data/rusers/sheddn/UCLA-ASD/data/ASD_BA9')


# ASD <- merge(BW1, y=c(BW3,BW4,BW7,BW9,BW10,BW14,BW16,BW17,BW18,BW19,BW21,BW22,BW24,BW26,BW27,BW28,BW29,BW31,BW32,BW34,
#                       BW37,BW38,BW39,BW42,BW44,BW45,BW48,BW49,BW54,BW59,BW60,BW65,BW66,BW67,BW68,BW70,BW72,BW73,BW75,BW76),
#              add.cell.ids=c('BW1','BW3','BW4','BW6','BW9','BW10','BW14','BW16','BW17','BW18','BW19','BW21','BW22','BW24',
#                             'BW26','BW27','BW28','BW29','BW31','BW32','BW34','BW37','BW38','BW39','BW42','BW44','BW45',
#                             'BW48','BW49','BW54','BW59','BW60','BW65','BW66','BW67','BW68','BW70','BW72','BW73','BW75','BW76'),
#              project = "UCLA-ASD")
# 
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_SampleLabels')
# 
# new.cluster.ids <- c('ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD',
#                      'ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD',
#                      'ASD','ASD','ASD','ASD','ASD','ASD','ASD')
# names(new.cluster.ids) <- levels(ASD)
# ASD <- RenameIdents(ASD, new.cluster.ids)
# 
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_GroupLabel')

