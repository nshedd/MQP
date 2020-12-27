library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#VisualCortex = readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

path1 = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

VisualCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

VisualCortex[["percent.mt"]] <- PercentageFeatureSet(VisualCortex, pattern = "^MT-")		

VisualCortex <- subset(VisualCortex, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

VisualCortex <- NormalizeData(VisualCortex, normalization.method = "LogNormalize", scale.factor = 10000)		

VisualCortex <- FindVariableFeatures(VisualCortex, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(VisualCortex)		
VisualCortex <- ScaleData(VisualCortex, features = all_cells)		

VisualCortex <- RunPCA(VisualCortex, features = VariableFeatures(object = VisualCortex), npcs=200)		

VisualCortex <- JackStraw(VisualCortex, dims=200, num.replicate = 100)		
VisualCortex <- ScoreJackStraw(VisualCortex, dims = 1:200)
Jackstraw_VisualCortex <- JackStrawPlot(VisualCortex, dims = 1:200)		
ggsave(path.expand("~/Lake/VisualCortex/jackstrawplot_GSE97930_VisualCortex_Seurat.png"), device=, width = 14, height = 7)		

Elbow_VisualCortex <- ElbowPlot(VisualCortex)		
ggsave(path.expand("~/Lake/VisualCortex/elbowplot_GSE97930_VisualCortex_Seurat.png"), device=, width = 14, height = 7)		

saveRDS(VisualCortex, file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))
