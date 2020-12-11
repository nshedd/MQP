library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

visualcortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

visualcortex <- NormalizeData(visualcortex, normalization.method = "LogNormalize", scale.factor = 10000)

visualcortex <- FindVariableFeatures(object = visualcortex)

all_cells <- rownames(visualcortex)
visualcortex <- ScaleData(visualcortex, features = all_cells)

visualcortex <- RunPCA(visualcortex, features = VariableFeatures(object = visualcortex))

visualcortex <- FindNeighbors(visualcortex, dims = 1:10)
visualcortex <- FindClusters(visualcortex, resolution = 0.5)

visualcortex <- RunUMAP(visualcortex, dims = 1:10)

new.cluster.ids <- c("Excitatory/Glutamatergic Neuron", "Neuron", "Oligodendrocyte", "GABAergic Neuron", "Glutamatergic Neuron", "Excitatory/Glutamatergic Neuron", 
                     "Astrocyte", "Non-myelinating Schwann Cell", "Neuron", "Oligodendrocyte Precursor", "GABAergic/Glutamatergic Neuron", "Microglia", "Endothelial/Pericyte")
names(new.cluster.ids) <- levels(visualcortex)
visualcortex <- RenameIdents(visualcortex, new.cluster.ids)

visualcortex.markers <- FindAllMarkers(visualcortex, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = visualcortex.markers %>% group_by(cluster)

saveRDS(visualcortex, file = path.expand("~/GSE97930_visualcortex_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(visualcortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/umap_GSE97930_visualcortex_Seurat_default.png"), device=)

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

write.table(diff_expressed, file = path.expand("~/GSE97930_visualcortex_differentiallyexpressed.txt"), sep="\t")
