library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


#CerebellarHem <- readRDS(file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

CerebellarHem <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

CerebellarHem[["percent.mt"]] <- PercentageFeatureSet(CerebellarHem, pattern = "^MT-")

CerebellarHem <- subset(CerebellarHem, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

CerebellarHem <- NormalizeData(CerebellarHem, normalization.method = "LogNormalize", scale.factor = 10000)

CerebellarHem <- FindVariableFeatures(object = CerebellarHem)

all_cells <- rownames(CerebellarHem)
CerebellarHem <- ScaleData(CerebellarHem, features = all_cells)

CerebellarHem <- RunPCA(CerebellarHem, features = VariableFeatures(object = CerebellarHem))

saveRDS(CerebellarHem, file = path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_Seurat.rds"))

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:20)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 0.5)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:20, metric="euclidean")

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

FeaturePlot(CerebellarHem, features = c("SLC17A7", "GRIN1", "GRIN2B"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_ExcitatoryFeatures.png"), device=)
 
FeaturePlot(CerebellarHem, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_InhibitoryFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_OligodendrocyteFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("GFAP", "SLC1A2", "SLC1A3", "SLC4A4", "AQP4"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_AstrocyteFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_EndothelialFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("COBLL1", "PDGFRB", "COBLL1", "PDGFRB"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_PericyteFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("PCDH15"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_OPCFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("APBB1IP", "P2RY12", "PTPRC"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_MicrogliaFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("RYR1"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_PurkinjeFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("RELN", "GRM4", "RBFOX3"))
ggsave(path.expand("~/Lake/CerebellarHem/GSE97930_CerebellarHem_Seurat_GranuleFeatures.png"), device=)
