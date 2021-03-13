library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)
library(reshape2)

##SingleR
## Load reference dataset
path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)
Lake <- CreateSeuratObject(counts = matrix, project = "SeuratPipeline", min.cells = 3, min.features = 200)

Lake[["percent.mt"]] <- PercentageFeatureSet(Lake, pattern = "^MT-")

Lake <- subset(Lake, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

Lake <- NormalizeData(Lake, normalization.method = "LogNormalize", scale.factor = 10000)

Lake_SCE <- as.SingleCellExperiment(Lake)
Lake_labels <- Idents(Lake)


## ASD-BA4.6
print("Loading ASD UMAP data w/o Doublets BA4.6...")
ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved_Relabeled.RDS')

ASD_SCE <- as.SingleCellExperiment(ASD)
ASD_clust <- Idents(ASD)

print('Running SingleR...')
ASD_SingleR <- SingleR(test=ASD_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=ASD_clust,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print('Plotting...')
ASD$SingleR.pruned.calls <- ASD_SingleR$pruned.labels

new_labels = c('NA')
for (i in ASD_SingleR$labels) {
  if (grepl("Ex", i, fixed=TRUE)){
    new_labels = c(new_labels, 'Ex')
  } else if (grepl("In", i, fixed=TRUE)){
    new_labels = c(new_labels, 'In')
  } else{
    new_labels = c(new_labels, i)
  }
}
ASD$SingleR.calls <- new_labels

ASD_BA4.6_prop_singler_list = ASD$SingleR.calls
ASD_BA4.6_prop_singler_table = table(ASD_BA4.6_prop_singler_list)

ASD_BA4.6_prop_df = ASD_BA4.6_prop_singler_table %>% as.data.frame
colnames(ASD_BA4.6_prop_df) <- c('Cell_Type', 'FreqASD')

ASD_BA4.6_sum = sum(ASD_BA4.6_prop_df$FreqASD)
ASD_BA4.6_prop_df$FreqASD = ASD_BA4.6_prop_df$FreqASD/ASD_BA4.6_sum


## CTL-BA4.6
print("Loading CTL UMAP data w/o Doublets BA4.6...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved_Relabeled.RDS')

CTL_SCE <- as.SingleCellExperiment(CTL)
CTL_clust <- Idents(CTL)

print('Running SingleR...')
CTL_SingleR <- SingleR(test=CTL_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=CTL_clust,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print('Plotting...')
CTL$SingleR.pruned.calls <- CTL_SingleR$pruned.labels

new_labels = c('NA')
for (i in CTL_SingleR$labels) {
  if (grepl("Ex", i, fixed=TRUE)){
    new_labels = c(new_labels, 'Ex')
  } else if (grepl("In", i, fixed=TRUE)){
    new_labels = c(new_labels, 'In')
  } else{
    new_labels = c(new_labels, i)
  }
}
CTL$SingleR.calls <- new_labels

CTL_BA4.6_prop_singler_list = CTL$SingleR.calls
CTL_BA4.6_prop_singler_table = table(CTL_BA4.6_prop_singler_list)

CTL_BA4.6_prop_df = CTL_BA4.6_prop_singler_table %>% as.data.frame
colnames(CTL_BA4.6_prop_df) <- c('Cell_Type', 'FreqCTL')

CTL_BA4.6_sum = sum(CTL_BA4.6_prop_df$FreqCTL)
CTL_BA4.6_prop_df$FreqCTL = CTL_BA4.6_prop_df$FreqCTL/CTL_BA4.6_sum


BA4.6_prop_df = merge(x = ASD_BA4.6_prop_df, y = CTL_BA4.6_prop_df, by = 'Cell_Type', all = TRUE)
print(BA4.6_prop_df)

write.table(BA4.6_prop_df, '/data/rusers/sheddn/UCLA-ASD/data/CellTypeProportions_SingleR_BA4.6.txt', sep=',')

BA4.6_prop_df = melt(BA4.6_prop_df, id.vars='Cell_Type')

ggplot(data=BA4.6_prop_df, aes(x = Cell_Type, y = value, fill=variable)) + geom_bar(stat="identity", position='dodge') +ggtitle("BA 4/6 - Cluster Labels")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CellTypeProportion_SingleR_BA4.6.png')


## ASD-BA9
print("Loading UMAP data w/o Doublets BA9...")
ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved_Relabeled.RDS')

ASD_SCE <- as.SingleCellExperiment(ASD)
ASD_clust <- Idents(ASD)

print('Running SingleR...')
ASD_SingleR <- SingleR(test=ASD_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=ASD_clust,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print('Plotting...')
ASD$SingleR.pruned.calls <- ASD_SingleR$pruned.labels

new_labels = c('NA')
for (i in ASD_SingleR$labels) {
  if (grepl("Ex", i, fixed=TRUE)){
    new_labels = c(new_labels, 'Ex')
  } else if (grepl("In", i, fixed=TRUE)){
    new_labels = c(new_labels, 'In')
  } else{
    new_labels = c(new_labels, i)
  }
}
ASD$SingleR.calls <- new_labels

ASD_BA9_prop_singler_list = ASD$SingleR.calls
ASD_BA9_prop_singler_table = table(ASD_BA9_prop_manual_list)

ASD_BA9_prop_df = ASD_BA9_prop_singler_table %>% as.data.frame
colnames(ASD_BA9_prop_df) <- c('Cell_Type', 'FreqASD')

ASD_BA9_sum = sum(ASD_BA9_prop_df$FreqASD)
ASD_BA9_prop_df$FreqASD = ASD_BA9_prop_df$FreqASD/ASD_BA9_sum


## CTL-BA9
print("Saving UMAP data w/o Doublets BA9...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved_Relabeled.RDS')

CTL_SCE <- as.SingleCellExperiment(CTL)
CTL_clust <- Idents(CTL)

print('Running SingleR...')
CTL_SingleR <- SingleR(test=CTL_SCE,
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=CTL_clust,
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")

print('Plotting...')
CTL$SingleR.pruned.calls <- CTL_SingleR$pruned.labels

new_labels = c('NA')
for (i in CTL_SingleR$labels) {
  if (grepl("Ex", i, fixed=TRUE)){
    new_labels = c(new_labels, 'Ex')
  } else if (grepl("In", i, fixed=TRUE)){
    new_labels = c(new_labels, 'In')
  } else{
    new_labels = c(new_labels, i)
  }
}
CTL$SingleR.calls <- new_labels

CTL_BA9_prop_singler_list = CTL$SingleR.calls
CTL_BA9_prop_singler_table = table(CTL_BA9_prop_singler_list)

CTL_BA9_prop_df = CTL_BA9_prop_singler_table %>% as.data.frame
colnames(CTL_BA9_prop_df) <- c('Cell_Type', 'FreqCTL')

CTL_BA9_sum = sum(CTL_BA9_prop_df$FreqCTL)
CTL_BA9_prop_df$FreqCTL = CTL_BA9_prop_df$FreqCTL/CTL_BA9_sum


BA9_prop_df = merge(x = ASD_BA9_prop_df, y = CTL_BA9_prop_df, by = 'Cell_Type', all = TRUE)
print(BA9_prop_df)

write.table(BA9_prop_df, '/data/rusers/sheddn/UCLA-ASD/data/CellTypeProportions_SingleR_BA9.txt', sep=',')

BA9_prop_df = melt(BA9_prop_df, id.vars='Cell_Type')

ggplot(data=BA9_prop_df, aes(x = Cell_Type, y = value, fill=variable)) + geom_bar(stat="identity", position='dodge') +ggtitle("BA 4/6 - Cluster Labels")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CellTypeProportion_SingleR_BA9.png')
