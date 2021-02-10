library(Seurat)
library(ggplot2)
library(dplyr)

path1 = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5376/matrix.tsv")
path2 = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5342/matrix.tsv")
path3 = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_UMB5168/matrix.tsv")

matrix1 = read.table(path1, header=TRUE, row.names=1, sep='\t')
matrix2 = read.table(path2, header=TRUE, row.names=1, sep='\t')
matrix3 = read.table(path3, header=TRUE, row.names=1, sep='\t')

combinedmatrix = merge(matrix1, matrix2, by=0)

row.names(combinedmatrix)<-combinedmatrix$Row.names
combinedmatrix$Row.names <- NULL

combinedmatrix = merge(combinedmatrix, matrix3, by=0)

row.names(combinedmatrix)<-combinedmatrix$Row.names
combinedmatrix$Row.names <- NULL

All <- CreateSeuratObject(counts = matrix, project = "PEC", min.cells = 3, min.features = 200)

All[["percent.mt"]] <- PercentageFeatureSet(All, pattern = "^MT-")

All <- NormalizeData(All, normalization.method = "LogNormalize", scale.factor = 10000)

All <- FindVariableFeatures(object = All)

all_cells <- rownames(All)
All <- ScaleData(All, features = all_cells)

All <- RunPCA(All, features = VariableFeatures(object = All))

ElbowPlot(All)
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/ElbowPlot.png"), device=)

#saveRDS(All, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/Matrix_Seurat.rds"))

All <- FindNeighbors(All, dims = 1:20)
All <- FindClusters(All, resolution = 0.5)

All <- RunUMAP(All, dims = 1:20, metric="euclidean", n.neighbors=10)

plot = DimPlot(All, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/umap_numberedclusters.png"), device=)

All.markers <- FindAllMarkers(All, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
All.markers %>% group_by(cluster)

brain_genes = read.table(path.expand("~/marker_genes.csv"), header=TRUE, sep=",")

diff_expressed = All.markers
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

write.table(diff_expressed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/diffexpressed_condensed.txt"), sep = '\t')

FeaturePlot(All, features = c("ENSG00000176884.15", "ENSG00000104888.10", "ENSG00000119042.17"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/ExcitatoryFeatures.png"), device=)

FeaturePlot(All, features = c("ENSG00000128683.14", "ENSG00000136750.13", "ENSG00000157103.12"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/InhibitoryFeatures.png"), device=)

FeaturePlot(All, features = c("ENSG00000204655.12", "ENSG00000013297.11", "ENSG00000168314.17"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/OligodendrocyteFeatures.png"), device=)

FeaturePlot(All, features = c("ENSG00000144908.13", "ENSG00000135821.18", "ENSG00000160307.10", "ENSG00000171885.17", "ENSG00000110436.13"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/AstrocyteFeatures.png"), device=)

FeaturePlot(All, features = c("ENSG00000082438.17"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/EndothelialFeatures.png"), device=)

FeaturePlot(All, features = c("ENSG00000134853.12", "ENSG00000184221.13"))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/OPCFeatures.png"), device=)

intersection = diff_expressed_condensed

dotplot <- DotPlot(All, features = intersection) + 
   theme(axis.text.x = element_text(angle = 90)) + 
   scale_y_discrete(limits = rev(levels(All$cluster)))
ggsave(path.expand("~/PEC_CTL_IsoHuB_DLPFC_snRNASeq_NextSeq500_All/analysis/dotplot.png"), device=)
