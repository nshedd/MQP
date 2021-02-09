library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)


path1 = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/matrix.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")

FrontalCortex <- subset(FrontalCortex, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)

FrontalCortex <- FindVariableFeatures(object = FrontalCortex)

all_cells <- rownames(FrontalCortex)
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))

saveRDS(FrontalCortex, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/Matrix_Seurat.rds"))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 1)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

plot = DimPlot(FrontalCortex, preduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/umap_numberedclusters.png"), device=)

FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FrontalCortex.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(brainRegion.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(brainRegion, features = intersection) + 
   theme(axis.text.x = element_text(angle = 90)) + 
   scale_y_discrete(limits = rev(levels(brainRegion$seurat_clusters)))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/dotplot.png"), device=)



brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = FrontalCortex.markers

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


write.table(diff_expressed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/diffexpressed_condensed.txt"), sep = '\t')

FeaturePlot(VisualCortex, features = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/ExcitatoryFeatures.png"), device=)

FeaturePlot(VisualCortex, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/InhibitoryFeatures.png"), device=)

FeaturePlot(VisualCortex, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/OligodendrocyteFeatures.png"), device=)

FeaturePlot(VisualCortex, features = c("SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/AstrocyteFeatures.png"), device=)

FeaturePlot(VisualCortex, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/EndothelialFeatures.png"), device=)

FeaturePlot(VisualCortex, features = c("PCDH15", "OLIG1", "PCDH15", "OLIG1"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/OPCFeatures.png"), device=)

FeaturePlot(VisualCortex, features = c("APBB1IP", "P2RY12", "APBB1IP", "P2RY12"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/MicrogliaFeatures.png"), device=)
