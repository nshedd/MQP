library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

#CerebellarHem = readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

CerebellarHem <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

CerebellarHem[["percent.mt"]] <- PercentageFeatureSet(CerebellarHem, pattern = "^MT-")		

CerebellarHem <- subset(CerebellarHem, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

CerebellarHem <- NormalizeData(CerebellarHem, normalization.method = "LogNormalize", scale.factor = 10000)		

CerebellarHem <- FindVariableFeatures(CerebellarHem, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(CerebellarHem)		
CerebellarHem <- ScaleData(CerebellarHem, features = all_cells)		

CerebellarHem <- RunPCA(CerebellarHem, features = VariableFeatures(object = CerebellarHem), npcs=200)		

CerebellarHem <- JackStraw(CerebellarHem, dims=200, num.replicate = 100)		
CerebellarHem <- ScoreJackStraw(CerebellarHem, dims = 1:200)
Jackstraw_CerebellarHem <- JackStrawPlot(CerebellarHem, dims = 1:200)		
ggsave(path.expand("~/Lake/CerebellarHem/jackstrawplot_GSE97930_CerebellarHem_Seurat.png"), device=, width = 14, height = 7)		

Elbow_CerebellarHem <- ElbowPlot(CerebellarHem)		
ggsave(path.expand("~/Lake/CerebellarHem/elbowplot_GSE97930_CerebellarHem_Seurat.png"), device=, width = 14, height = 7)		

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:2)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:2, metric="euclidean")

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(CerebellarHem)
#CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

CerebellarHem.markers <- FindAllMarkers(CerebellarHem, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = CerebellarHem.markers %>% group_by(cluster)

saveRDS(CerebellarHem, file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_02pc.png"), device=)
