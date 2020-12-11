library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

frontalcortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

frontalcortex <- NormalizeData(frontalcortex, normalization.method = "LogNormalize", scale.factor = 10000)

frontalcortex <- FindVariableFeatures(object = frontalcortex)

all_cells <- rownames(frontalcortex)
frontalcortex <- ScaleData(frontalcortex, features = all_cells)

frontalcortex <- RunPCA(frontalcortex, features = VariableFeatures(object = frontalcortex))

frontalcortex <- FindNeighbors(frontalcortex, dims = 1:10)
frontalcortex <- FindClusters(frontalcortex, resolution = 0.5)

frontalcortex <- RunUMAP(frontalcortex, dims = 1:10)

new.cluster.ids <- c("Oligodendrocyte", "Neuron", "Excitatory Neuron", "GABAergic Neuron", "Glutamatergic Neuron", "GABAergic/Glutamatergic Neuron", 
    "Astrocyte", "NA", "Neuron", "Oligodendrocyte Precursor", "Neuron", "Microglia")
names(new.cluster.ids) <- levels(frontalcortex)
frontalcortex <- RenameIdents(frontalcortex, new.cluster.ids)

frontalcortex.markers <- FindAllMarkers(frontalcortex, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = frontalcortex.markers %>% group_by(cluster)

saveRDS(frontalcortex, file = path.expand("~/GSE97930_frontalcortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(frontalcortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/umap_GSE97930_frontalcortex_Seurat_default.png"), device=)

path2 = path.expand("~/Zlab single-cell marker genes - Brain.tsv")

brain_genes = read.table(path2, header=TRUE, sep= "\t")
blood_genes = 
  
  celltypes <- character()
for (gene in diff_expressed$gene) {
  if (gene %in% brain_genes$Human_Gene) {
    celltypes <- c(celltypes, T_cells$Human_Gene)
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

write.table(diff_expressed, file = path.expand("~/GSE97930_frontalcortex_differentiallyexpressed.txt"), sep="\t")
