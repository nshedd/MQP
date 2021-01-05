library(Seurat)
library(ggplot2)

FrontalCortex <- readRDS(file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:35)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:35, metric="euclidean")

#new.cluster.ids <- c("In1","Ex1","Ex2","End","Oli1","Ex3","Oli1","Ex4","Ex5","Per/Mic","OPC","Ex6","Ex7","In2","In3","Ast","Oli2")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_035pc.png"), device=)
