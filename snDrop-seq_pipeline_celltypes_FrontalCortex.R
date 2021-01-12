library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)		

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")		

FrontalCortex <- subset(FrontalCortex, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)		

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)		

FrontalCortex <- FindVariableFeatures(FrontalCortex, selection.method = "vst", nfeatures = 2000)		

all_cells <- rownames(FrontalCortex)		
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)		

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))		

saveRDS(FrontalCortex, file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 1.0)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_findct.png"), device=)

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster)

brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

head(brain_genes)

diff_expressed = pbmc.markers

head(diff_expressed)

celltypelist = brain_genes$Cell.type

celltypes <- character()
for (gene in diff_expressed$name) {
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


write.table(diff_expressed, file = path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_diffexpressed_condensed.txt"), sep = '\t')

#FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))
