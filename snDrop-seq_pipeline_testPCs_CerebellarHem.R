library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

CerebellarHem = readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:5)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:5, metric="euclidean")

new.cluster.ids <- c("Gran1", "Gran2", "Gran3", "Purk1", "Ast", "Purk2", "OPC", "Mic/End/Per", "Oli")
names(new.cluster.ids) <- levels(CerebellarHem)
CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

CerebellarHem.markers <- FindAllMarkers(CerebellarHem, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = CerebellarHem.markers %>% group_by(cluster)


plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_005pc.png"), device=)
