devtools::install_github("nshedd/scclusteval")

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scclusteval)



## Paper labels

# Load RDS from CellTypeIdentification Script
All2 <- readRDS(file=path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

# Run and plot UMAP without generatingclusters
All2 <- RunUMAP(All2, dims = 1:20, metric="euclidean")

plot = DimPlot(All2, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_findct_oglabels.png"), device=)

## My labels

# Load RDS from CellTypeIdentification Script
All <- readRDS(file=path.expand("~/Lake/All/GSE97930_All_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

# Find Neighbors and CLusters again (it will generate the same clusters)
All <- FindNeighbors(All, dims = 1:20)
All <- FindClusters(All, resolution = 1.5)

All <- RunUMAP(All, dims = 1:20, metric="euclidean")

ggsave(path.expand("~/Lake/All/umap_GSE97930_All_Seurat_findct.png"), device=)

# Manually relabel clusters based on analysis
new.cluster.ids <- c("Neuron", "Oli1", "Ex1", "Ex2", "In1", "Ast1", "Ex3", "Ex4", "In2", "In3",
                    "OPC", "Ex5", "Gran1", "Ex6", "Ex7", "Oli2", "Gran2", "Purk1", "Ex8", "Gran3",
                    "In4", "In5", "Mic", "Ex9", "Ast2", "Ex10", "End/Per", "Ex11", "Ast/Oli", "Mic/Oli", "Ex12", "Oli3")
names(new.cluster.ids) <- levels(All)
All <- RenameIdents(All, new.cluster.ids)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()



## Heatmap comparing labels

# Create factor vectors from idents to use in heatmap
oglabels <- Idents(All2)
oglabels <- factor(oglabels)

newlabels <- Idents(All)
newlabels <- factor(newlabels)

# Called from scclusteval script that I took off github
png(file=path.expand("~/Lake/All/GSE97930_All_Seurat_JaccardHeatmap.png"))
PairWiseJaccardSetsHeatmap(ident1 = oglabels, ident2 = newlabels)
dev.off()
