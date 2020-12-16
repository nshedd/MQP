library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


matrix = readRDS(path.expand("~/GSE97930_All.RDS"))

lake_all <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

lake_all <- NormalizeData(lake_all, normalization.method = "LogNormalize", scale.factor = 10000)

lake_all <- FindVariableFeatures(object = lake_all)

all_cells <- rownames(lake_all)
lake_all <- ScaleData(lake_all, features = all_cells)

lake_all <- RunPCA(lake_all, features = VariableFeatures(object = lake_all))

lake_all <- FindNeighbors(lake_all, dims = 1:10)
lake_all <- FindClusters(lake_all, resolution = 0.5)

new.cluster.ids <- c("Neuron 1", "Neuron 2", "Oligodendrocyte", "Excitatory Neuron 1", "Inhibitory Neuron 1", "Inhibitory Neuron 2", "Neuron 3",
                     "Astrocyte 1", "Excitatory Neuron 2", "Oligodendrocyte Precursor", "Inhibitory Neuron 3", "? 1", "Microglia", "Astrocyte 2",
                     "Endothelial", "? 2", "Astrocyte 3")

names(new.cluster.ids) <- levels(lake_all)
lake_all <- RenameIdents(lake_all, new.cluster.ids)

lake_all <- RunUMAP(lake_all, dims = 1:10)

saveRDS(lake_all, file = path.expand("~/GSE97930_lake_all_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(lake_all, reduction = "umap")
ggsave(path.expand("~/umap_GSE97930_lake_all_Seurat_default.png"), device=)

lake_all.markers <- FindAllMarkers(lake_all, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = lake_all.markers %>% group_by(cluster)

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

write.table(diff_expressed, file = path.expand("~/GSE97930_lake_all_differentiallyexpressed.txt"), sep="\t")
