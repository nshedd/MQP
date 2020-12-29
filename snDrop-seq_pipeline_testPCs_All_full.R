library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


matrix = matrix = readRDS(path.expand("~/GSE97930_All.RDS"))

All <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

All[["percent.mt"]] <- PercentageFeatureSet(All, pattern = "^MT-")		

All <- subset(All, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

All <- NormalizeData(All, normalization.method = "LogNormalize", scale.factor = 10000)		

All <- FindVariableFeatures(All, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(All)		
All <- ScaleData(All, features = all_cells)		

All <- RunPCA(All, features = VariableFeatures(object = All), npcs=200)		

All <- JackStraw(All, dims=200, num.replicate = 100)		
All <- ScoreJackStraw(All, dims = 1:200)
Jackstraw_All <- JackStrawPlot(All, dims = 1:200)		
ggsave(path.expand("~/Lake/All/jackstrawplot_GSE97930_All_Seurat.png"), device=, width = 14, height = 7)		

Elbow_All <- ElbowPlot(All)		
ggsave(path.expand("~/Lake/All/elbowplot_GSE97930_All_Seurat.png"), device=, width = 14, height = 7)		

saveRDS(All, file = path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))
