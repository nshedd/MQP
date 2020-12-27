library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

CerebellarHem = readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:2, metric="euclidean")

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(CerebellarHem)
#CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

CerebellarHem.markers <- FindAllMarkers(CerebellarHem, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = CerebellarHem.markers %>% group_by(cluster)

saveRDS(CerebellarHem, file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_02pc_short.png"), device=)
