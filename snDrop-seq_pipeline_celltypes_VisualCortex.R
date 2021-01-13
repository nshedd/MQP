library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


#VisualCortex <- readRDS(file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

path1 = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

VisualCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

VisualCortex[["percent.mt"]] <- PercentageFeatureSet(VisualCortex, pattern = "^MT-")

VisualCortex <- NormalizeData(VisualCortex, normalization.method = "LogNormalize", scale.factor = 10000)

VisualCortex <- FindVariableFeatures(object = VisualCortex)

all_cells <- rownames(VisualCortex)
VisualCortex <- ScaleData(VisualCortex, features = all_cells)

VisualCortex <- RunPCA(VisualCortex, features = VariableFeatures(object = VisualCortex))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:10)
VisualCortex <- FindClusters(VisualCortex, resolution = 0.5)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:10)

saveRDS(VisualCortex, file = path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(VisualCortex)
#VisualCortex <- RenameIdents(VisualCortex, new.cluster.ids)

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_findct.png"), device=)

VisualCortex.markers <- FindAllMarkers(VisualCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VisualCortex.markers %>% group_by(cluster)

brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = VisualCortex.markers

celltypelist = brain_genes$Cell.type

celltypes <- character()
for (gene in diff_expressed$gene) {
  if (gene %in% brain_genes$Human.Gene) {
    celltype = brain_genes$Cell.type[brain_genes$Human.Gene == gene]
    celltype = paste(celltype, collapse = ', ')
    celltypes <- c(celltypes, celltype)
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

diff_expressed_condensed = diff_expressed[diff_expressed$cell_type != "unknown",]


write.table(diff_expressed, file = path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/Lake/VisualCortex/umap_GSE97930_VisualCortex_Seurat_diffexpressed_condensed.txt"), sep = '\t')

#FeaturePlot(VisualCortex, features = c("APBB1IP", "P2RY12"))
#ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_MicrogliaFeatures.png"), device=)
