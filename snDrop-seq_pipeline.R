library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

CerebellarHem <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

CerebellarHem <- NormalizeData(CerebellarHem, normalization.method = "LogNormalize", scale.factor = 10000)

CerebellarHem <- FindVariableFeatures(object = CerebellarHem)

all_cells <- rownames(CerebellarHem)
CerebellarHem <- ScaleData(CerebellarHem, features = all_cells)

CerebellarHem <- RunPCA(CerebellarHem, features = VariableFeatures(object = CerebellarHem))

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:10)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:10)

new.cluster.ids <- c("Neuron", "Neuron", "Neuron", "Neuron", "Astrocyte", "Neuron", "?", "Oligodendrocyte", "Microglia", "Neuron", "?", "?")
names(new.cluster.ids) <- levels(CerebellarHem)
CerebellarHem <- RenameIdents(CerebellarHem, new.cluster.ids)

CerebellarHem.markers <- FindAllMarkers(CerebellarHem, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = CerebellarHem.markers %>% group_by(cluster)

saveRDS(CerebellarHem, file = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/umap_GSE97930_CerebellarHem_Seurat_default.png"), device=)

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

write.table(diff_expressed, file = path.expand("~/GSE97930_CerebellarHem_differentiallyexpressed.txt"), sep="\t")
