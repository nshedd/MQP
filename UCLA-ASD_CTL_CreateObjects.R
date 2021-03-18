library(dplyr)
library(Seurat)
library(ggplot2)

print('2BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/2BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW2 = CreateSeuratObject(counts = expression_matrix, project='2BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
 if(min(BW2$nCount_RNA)>=minCov){
    countLOW=min(BW2$nCount_RNA)
  }else{
    countLOW=quantile(BW2$nCount_RNA, prob=c(0.01))  
 }
 countHIGH=quantile(BW2$nCount_RNA, prob=0.97) 
 featureLOW=quantile(BW2$nFeature_RNA, prob=0.01)
BW2 <- subset(BW2, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('5BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/5BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW5 = CreateSeuratObject(counts = expression_matrix, project='5BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW5$nCount_RNA)>=minCov){
  countLOW=min(BW5$nCount_RNA)
}else{
  countLOW=quantile(BW5$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW5$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW5$nFeature_RNA, prob=0.01)
BW5 <- subset(BW5, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('6BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/6BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW6 = CreateSeuratObject(counts = expression_matrix, project='6BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW6$nCount_RNA)>=minCov){
  countLOW=min(BW6$nCount_RNA)
}else{
  countLOW=quantile(BW6$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW6$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW6$nFeature_RNA, prob=0.01)
BW6 <- subset(BW6, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('8BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/8BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW8 = CreateSeuratObject(counts = expression_matrix, project='8BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW8$nCount_RNA)>=minCov){
  countLOW=min(BW8$nCount_RNA)
}else{
  countLOW=quantile(BW8$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW8$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW8$nFeature_RNA, prob=0.01)
BW8 <- subset(BW8, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('11BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/11BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW11 = CreateSeuratObject(counts = expression_matrix, project='11BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW11$nCount_RNA)>=minCov){
  countLOW=min(BW11$nCount_RNA)
}else{
  countLOW=quantile(BW11$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW11$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW11$nFeature_RNA, prob=0.01)
BW11 <- subset(BW11, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('12BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/12BW_S15/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW12 = CreateSeuratObject(counts = expression_matrix, project='12BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW12$nCount_RNA)>=minCov){
  countLOW=min(BW12$nCount_RNA)
}else{
  countLOW=quantile(BW12$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW12$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW12$nFeature_RNA, prob=0.01)
BW12 <- subset(BW12, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('13BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/13BW_S16/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW13 = CreateSeuratObject(counts = expression_matrix, project='13BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW13$nCount_RNA)>=minCov){
  countLOW=min(BW13$nCount_RNA)
}else{
  countLOW=quantile(BW13$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW13$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW13$nFeature_RNA, prob=0.01)
BW13 <- subset(BW13, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('15BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/15BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW15 = CreateSeuratObject(counts = expression_matrix, project='15BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW15$nCount_RNA)>=minCov){
  countLOW=min(BW15$nCount_RNA)
}else{
  countLOW=quantile(BW15$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW15$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW15$nFeature_RNA, prob=0.01)
BW15 <- subset(BW15, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('20BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/20BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW20 = CreateSeuratObject(counts = expression_matrix, project='20BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW20$nCount_RNA)>=minCov){
  countLOW=min(BW20$nCount_RNA)
}else{
  countLOW=quantile(BW20$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW20$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW20$nFeature_RNA, prob=0.01)
BW20 <- subset(BW20, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('23BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/23BW_S9/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW23 = CreateSeuratObject(counts = expression_matrix, project='23BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW23$nCount_RNA)>=minCov){
  countLOW=min(BW23$nCount_RNA)
}else{
  countLOW=quantile(BW23$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW23$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW23$nFeature_RNA, prob=0.01)
BW23 <- subset(BW23, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('25BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/25BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW25 = CreateSeuratObject(counts = expression_matrix, project='25BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW25$nCount_RNA)>=minCov){
  countLOW=min(BW25$nCount_RNA)
}else{
  countLOW=quantile(BW25$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW25$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW25$nFeature_RNA, prob=0.01)
BW25 <- subset(BW25, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('30BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/30BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW30 = CreateSeuratObject(counts = expression_matrix, project='30BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW30$nCount_RNA)>=minCov){
  countLOW=min(BW30$nCount_RNA)
}else{
  countLOW=quantile(BW30$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW30$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW30$nFeature_RNA, prob=0.01)
BW30 <- subset(BW30, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('33BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/33BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW33 = CreateSeuratObject(counts = expression_matrix, project='33BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW33$nCount_RNA)>=minCov){
  countLOW=min(BW33$nCount_RNA)
}else{
  countLOW=quantile(BW33$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW33$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW33$nFeature_RNA, prob=0.01)
BW33 <- subset(BW33, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('35BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/35BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW35 = CreateSeuratObject(counts = expression_matrix, project='35BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW35$nCount_RNA)>=minCov){
  countLOW=min(BW35$nCount_RNA)
}else{
  countLOW=quantile(BW35$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW35$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW35$nFeature_RNA, prob=0.01)
BW35 <- subset(BW35, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('40BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/40BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW40 = CreateSeuratObject(counts = expression_matrix, project='40BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW40$nCount_RNA)>=minCov){
  countLOW=min(BW40$nCount_RNA)
}else{
  countLOW=quantile(BW40$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW40$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW40$nFeature_RNA, prob=0.01)
BW40 <- subset(BW40, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('41BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/41BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW41 = CreateSeuratObject(counts = expression_matrix, project='41BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW41$nCount_RNA)>=minCov){
  countLOW=min(BW41$nCount_RNA)
}else{
  countLOW=quantile(BW41$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW41$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW41$nFeature_RNA, prob=0.01)
BW41 <- subset(BW41, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('43BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/43BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW43 = CreateSeuratObject(counts = expression_matrix, project='43BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW43$nCount_RNA)>=minCov){
  countLOW=min(BW43$nCount_RNA)
}else{
  countLOW=quantile(BW43$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW43$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW43$nFeature_RNA, prob=0.01)
BW43 <- subset(BW43, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('46BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/46BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW46 = CreateSeuratObject(counts = expression_matrix, project='46BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW46$nCount_RNA)>=minCov){
  countLOW=min(BW46$nCount_RNA)
}else{
  countLOW=quantile(BW46$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW46$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW46$nFeature_RNA, prob=0.01)
BW46 <- subset(BW46, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('47BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/47BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW47 = CreateSeuratObject(counts = expression_matrix, project='47BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW47$nCount_RNA)>=minCov){
  countLOW=min(BW47$nCount_RNA)
}else{
  countLOW=quantile(BW47$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW47$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW47$nFeature_RNA, prob=0.01)
BW47 <- subset(BW47, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('51BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/51BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW51 = CreateSeuratObject(counts = expression_matrix, project='51BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW51$nCount_RNA)>=minCov){
  countLOW=min(BW51$nCount_RNA)
}else{
  countLOW=quantile(BW51$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW51$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW51$nFeature_RNA, prob=0.01)
BW51 <- subset(BW51, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('52BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/52BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW52 = CreateSeuratObject(counts = expression_matrix, project='52BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW52$nCount_RNA)>=minCov){
  countLOW=min(BW52$nCount_RNA)
}else{
  countLOW=quantile(BW52$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW52$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW52$nFeature_RNA, prob=0.01)
BW52 <- subset(BW52, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('53BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/53BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW53 = CreateSeuratObject(counts = expression_matrix, project='53BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW53$nCount_RNA)>=minCov){
  countLOW=min(BW53$nCount_RNA)
}else{
  countLOW=quantile(BW53$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW53$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW53$nFeature_RNA, prob=0.01)
BW53 <- subset(BW53, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('55BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/55BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW55 = CreateSeuratObject(counts = expression_matrix, project='55BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW55$nCount_RNA)>=minCov){
  countLOW=min(BW55$nCount_RNA)
}else{
  countLOW=quantile(BW55$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW55$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW55$nFeature_RNA, prob=0.01)
BW55 <- subset(BW55, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('56BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/56BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW56 = CreateSeuratObject(counts = expression_matrix, project='56BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW56$nCount_RNA)>=minCov){
  countLOW=min(BW56$nCount_RNA)
}else{
  countLOW=quantile(BW56$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW56$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW56$nFeature_RNA, prob=0.01)
BW56 <- subset(BW56, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('57BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/57BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW57 = CreateSeuratObject(counts = expression_matrix, project='57BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW57$nCount_RNA)>=minCov){
  countLOW=min(BW57$nCount_RNA)
}else{
  countLOW=quantile(BW57$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW57$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW57$nFeature_RNA, prob=0.01)
BW57 <- subset(BW57, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('58BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/58BW_S28/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW58 = CreateSeuratObject(counts = expression_matrix, project='58BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW58$nCount_RNA)>=minCov){
  countLOW=min(BW58$nCount_RNA)
}else{
  countLOW=quantile(BW58$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW58$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW58$nFeature_RNA, prob=0.01)
BW58 <- subset(BW58, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('61BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/61BW_S31/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW61 = CreateSeuratObject(counts = expression_matrix, project='61BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW61$nCount_RNA)>=minCov){
  countLOW=min(BW61$nCount_RNA)
}else{
  countLOW=quantile(BW61$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW61$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW61$nFeature_RNA, prob=0.01)
BW61 <- subset(BW61, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('62BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/62BW_S32/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW62 = CreateSeuratObject(counts = expression_matrix, project='62BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW62$nCount_RNA)>=minCov){
  countLOW=min(BW62$nCount_RNA)
}else{
  countLOW=quantile(BW62$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW62$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW62$nFeature_RNA, prob=0.01)
BW62 <- subset(BW62, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('63BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/63BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW63 = CreateSeuratObject(counts = expression_matrix, project='63BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW63$nCount_RNA)>=minCov){
  countLOW=min(BW63$nCount_RNA)
}else{
  countLOW=quantile(BW63$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW63$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW63$nFeature_RNA, prob=0.01)
BW63 <- subset(BW63, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('64BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/64BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW64 = CreateSeuratObject(counts = expression_matrix, project='64BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW64$nCount_RNA)>=minCov){
  countLOW=min(BW64$nCount_RNA)
}else{
  countLOW=quantile(BW64$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW64$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW64$nFeature_RNA, prob=0.01)
BW64 <- subset(BW64, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('69BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/69BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW69 = CreateSeuratObject(counts = expression_matrix, project='69BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW69$nCount_RNA)>=minCov){
  countLOW=min(BW69$nCount_RNA)
}else{
  countLOW=quantile(BW69$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW69$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW69$nFeature_RNA, prob=0.01)
BW69 <- subset(BW69, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('71BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/71BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW71 = CreateSeuratObject(counts = expression_matrix, project='71BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW71$nCount_RNA)>=minCov){
  countLOW=min(BW71$nCount_RNA)
}else{
  countLOW=quantile(BW71$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW71$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW71$nFeature_RNA, prob=0.01)
BW71 <- subset(BW71, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('74BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/74BW_S16/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW74 = CreateSeuratObject(counts = expression_matrix, project='74BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW74$nCount_RNA)>=minCov){
  countLOW=min(BW74$nCount_RNA)
}else{
  countLOW=quantile(BW74$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW74$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW74$nFeature_RNA, prob=0.01)
BW74 <- subset(BW74, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('77BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/77BW_S19/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW77 = CreateSeuratObject(counts = expression_matrix, project='77BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW77$nCount_RNA)>=minCov){
  countLOW=min(BW77$nCount_RNA)
}else{
  countLOW=quantile(BW77$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW77$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW77$nFeature_RNA, prob=0.01)
BW77 <- subset(BW77, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

print('78BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/78BW_S20/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW78 = CreateSeuratObject(counts = expression_matrix, project='78BW')
minCov=1000 #if a sample has a good coverage (>=minCov), then don't set a lower thresold for nCount, it's already pretty good. 
if(min(BW78$nCount_RNA)>=minCov){
  countLOW=min(BW78$nCount_RNA)
}else{
  countLOW=quantile(BW78$nCount_RNA, prob=c(0.01))  
}
countHIGH=quantile(BW78$nCount_RNA, prob=0.97) 
featureLOW=quantile(BW78$nFeature_RNA, prob=0.01)
BW78 <- subset(BW78, subset = nFeature_RNA > featureLOW & nCount_RNA > countLOW  & nCount_RNA < countHIGH)

CTL_BA4.6 <- merge(BW11, y=c(BW12,BW13,BW15,BW2,BW23,BW25,BW30,BW35,BW46,BW5,BW6,BW8),
                    add.cell.ids=c('BW11','BW12','BW13','BW15','BW2','BW23','BW25','BW30','BW35','BW46','BW5','BW6','BW8'),
                   project = "UCLA-ASD")

saveRDS(CTL_BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6')

CTL_BA9 <- merge(BW20, y=c(BW33,BW40,BW41,BW43,BW47,BW51,BW52,BW53,BW55,BW56,BW57,BW61,BW62,BW63,BW64,BW69,BW71,BW74,BW77,BW78),
                    add.cell.ids=c('BW20','BW33','BW40','BW41','BW43','BW47','BW51','BW52','BW53','BW55','BW56','BW57','BW61','BW62','BW63','BW64','BW69','BW71','BW74','BW77','BW78'),
                   project = "UCLA-ASD")

saveRDS(CTL_BA9, '/data/rusers/sheddn/UCLA-ASD/data/CTL_BA9')


# CTL <- merge(BW2, y=c(BW5,BW6,BW8,BW11,BW12,BW13,BW15,BW20,BW23,BW25,BW30,BW33,BW35,BW40,BW41,BW43,BW46,BW47,
#                       BW51,BW52,BW53,BW55,BW56,BW57,BW58,BW61,BW62,BW63,BW64,BW69,BW71,BW74,BW77,BW78),
#              add.cell.ids=c('BW2','BW5','BW6','BW8','BW11','BW12','BW13','BW15','BW20','BW23','BW25','BW30','BW33','BW35','BW40','BW41','BW43','BW46',
#                             'BW47','BW51','BW52','BW53','BW55','BW56','BW57','BW58','BW61','BW62','BW63','BW64','BW69','BW71','BW74','BW77','BW78'),
#              project = "UCLA-ASD")
# 
# 
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_SampleLabels')
# 
# new.cluster.ids <- c('CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL',
#                      'CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL', 'CTL')
# names(new.cluster.ids) <- levels(CTL)
# CTL <- RenameIdents(CTL, new.cluster.ids)
# 
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_GroupLabel')


