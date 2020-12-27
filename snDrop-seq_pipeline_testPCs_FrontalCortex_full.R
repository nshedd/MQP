library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#FrontalCortex = readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")		

FrontalCortex <- subset(FrontalCortex, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)		

FrontalCortex <- FindVariableFeatures(FrontalCortex, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(FrontalCortex)		
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)		

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex), npcs=200)		

FrontalCortex <- JackStraw(FrontalCortex, dims=200, num.replicate = 100)		
FrontalCortex <- ScoreJackStraw(FrontalCortex, dims = 1:200)
Jackstraw_FrontalCortex <- JackStrawPlot(FrontalCortex, dims = 1:200)		
ggsave(path.expand("~/Lake/FrontalCortex/jackstrawplot_GSE97930_FrontalCortex_Seurat.png"), device=, width = 14, height = 7)		

Elbow_FrontalCortex <- ElbowPlot(FrontalCortex)		
ggsave(path.expand("~/Lake/FrontalCortex/elbowplot_GSE97930_FrontalCortex_Seurat.png"), device=, width = 14, height = 7)		

saveRDS(FrontalCortex, file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:2)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:2, metric="euclidean")

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_002pc.png"), device=)
