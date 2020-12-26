library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

CerebellarHem = readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:200)
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

#featureplot_CerebellarHem < - FeaturePlot(CerebellarHem, features = c("SYT1", "SYT1", "RBFOX3", "GAD2", "SLC1A3", "GRIN2B", "PCDH15", "MBP", "APBB1IP", "SLC6A1", "PCDH15", "SLC1A3"))
#ggsave(path.expand("~/Lake/CerebellarHem/featureplot_GSE97930_CerebellarHem_Seurat_default.png"), device=)
  
#write.table(diff_expressed, file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_differentiallyexpressed_origident.txt"), sep="\t")
