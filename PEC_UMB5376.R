library(Seurat)
library(ggplot2)
library(dplyr)

print("loaded libraries")

path1 = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/matrix.tsv")

matrix = read.table(path1, header=TRUE, row.names=1, sep='\t')

UMB5376 <- CreateSeuratObject(counts = matrix, project = "PEC", min.cells = 3, min.features = 200)

UMB5376[["percent.mt"]] <- PercentageFeatureSet(UMB5376, pattern = "^MT-")

UMB5376 <- NormalizeData(UMB5376, normalization.method = "LogNormalize", scale.factor = 10000)

UMB5376 <- FindVariableFeatures(object = UMB5376)

all_cells <- rownames(UMB5376)
UMB5376 <- ScaleData(UMB5376, features = all_cells)

UMB5376 <- RunPCA(UMB5376, features = VariableFeatures(object = UMB5376))

ElbowPlot(UMB5376)
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/ElbowPlot.png"), device=)

#saveRDS(UMB5376, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/Matrix_Seurat.rds"))

UMB5376 <- FindNeighbors(UMB5376, dims = 1:20)
UMB5376 <- FindClusters(UMB5376, resolution = 0.5)

UMB5376 <- RunUMAP(UMB5376, dims = 1:20, metric="euclidean", n.neighbors=10)

plot = DimPlot(UMB5376, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/umap_numberedclusters.png"), device=)

UMB5376.markers <- FindAllMarkers(UMB5376, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
UMB5376.markers %>% group_by(cluster)

brain_genes = read.table(path.expand("~/marker_genes.csv"), header=TRUE, sep=",")

diff_expressed = UMB5376.markers
diff_expressed_genes_short = substr(diff_expressed$gene,1,nchar(diff_expressed$gene)-3) 

celltypelist = brain_genes$Cell.type

celltypes <- character()
for (gene in diff_expressed_genes_short) {
   if (gene %in% brain_genes$Full.Gene.Name) {
      celltype = brain_genes$Cell.type[brain_genes$Full.Gene.Name == gene]
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

FeaturePlot(UMB5376, features = c("ENSG00000176884.15", "ENSG00000104888.10", "ENSG00000119042.17"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/ExcitatoryFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("ENSG00000128683.14", "ENSG00000136750.13", "ENSG00000157103.12"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/InhibitoryFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("ENSG00000204655.12", "ENSG00000013297.11", "ENSG00000168314.17"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/OligodendrocyteFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("ENSG00000144908.13", "ENSG00000135821.18", "ENSG00000160307.10", "ENSG00000171885.17", "ENSG00000110436.13"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/AstrocyteFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("ENSG00000082438.17"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/EndothelialFeatures.png"), device=)

FeaturePlot(UMB5376, features = c("ENSG00000134853.12", "ENSG00000184221.13"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/OPCFeatures.png"), device=)

intersection = diff_expressed_condensed

dotplot <- DotPlot(UMB5376, features = intersection) + 
   theme(axis.text.x = element_text(angle = 90)) + 
   scale_y_discrete(limits = rev(levels(UMB5376$seurat_clusters)))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/analysis/dotplot.png"), device=)
