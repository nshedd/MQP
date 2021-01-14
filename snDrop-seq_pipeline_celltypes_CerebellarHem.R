library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


#CerebellarHem <- readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

CerebellarHem <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

CerebellarHem[["percent.mt"]] <- PercentageFeatureSet(CerebellarHem, pattern = "^MT-")

CerebellarHem <- NormalizeData(CerebellarHem, normalization.method = "LogNormalize", scale.factor = 10000)

CerebellarHem <- FindVariableFeatures(object = CerebellarHem)

all_cells <- rownames(CerebellarHem)
CerebellarHem <- ScaleData(CerebellarHem, features = all_cells)

CerebellarHem <- RunPCA(CerebellarHem, features = VariableFeatures(object = CerebellarHem))

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:10)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:10, metric="euclidean")

saveRDS(CerebellarHem, file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))


plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_findct.png"), device=)

CerebellarHem.markers <- FindAllMarkers(CerebellarHem, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CerebellarHem.markers %>% group_by(cluster)

brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = CerebellarHem.markers

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


write.table(diff_expressed, file = path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/Lake/CerebellarHem/umap_GSE97930_CerebellarHem_Seurat_diffexpressed_condensed.txt"), sep = '\t')

# FeaturePlot(VisualCortex, features = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B"))
# ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_ExcitatoryFeatures.png"), device=)
# 
# FeaturePlot(VisualCortex, features = c("GAD1", "GAD2", "SLC6A1"))
# ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_InhibitoryFeatures.png"), device=)
# 
# FeaturePlot(VisualCortex, features = c("CLDN11", "MOG", "MOBP", "MBP"))
# ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_OligodendrocyteFeatures.png"), device=)
# 
# FeaturePlot(VisualCortex, features = c("SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4"))
# ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_AstrocyteFeatures.png"), device=)
# 
# FeaturePlot(VisualCortex, features = c("COBLL1", "DUSP1", "FLT1"))
# ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_EndothelialFeatures.png"), device=)
# 
# FeaturePlot(VisualCortex, features = c("PCDH15", "OLIG1", "PCDH15", "OLIG1"))
# ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_OPCFeatures.png"), device=)
# 
# FeaturePlot(VisualCortex, features = c("APBB1IP", "P2RY12", "APBB1IP", "P2RY12"))
# ggsave(path.expand("~/Lake/VisualCortex/GSE97930_VisualCortex_Seurat_MicrogliaFeatures.png"), device=)
