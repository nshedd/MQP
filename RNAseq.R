library(Seurat)
library(ggplot2)

VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:44)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:44, metric="euclidean")

new.cluster.ids <- c("Ex2","Ex1/In1","Oli","Ex3","In2","Ex4","Ast","In3","Ex5","In4","OPC","In5","Ex6","Mic","Ex7","End/Per","In6","Ex8")
names(new.cluster.ids) <- levels(VisualCortex)
VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_044pc.png"), device=)
