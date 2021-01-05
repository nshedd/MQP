library(Seurat)
library(ggplot2)

FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:35)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:35, metric="euclidean")

new.cluster.ids <- c("End/Ex1/In1","Ex2","Oli","Ex3","In2","Ast","In3","Ex4","In4","OPC","Ex5","In5","Mic","Ex6","Ex7","In6","Ex8","Per")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_035pc_oglabels.png"), device=)
