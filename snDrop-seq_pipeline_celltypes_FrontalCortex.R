library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DoubletFinder)


path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "set2", min.cells = 3, min.features = 200)

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")

FrontalCortex <- subset(FrontalCortex, subset = nFeature_RNA > 200 & nFeature_RNA < 3000)

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)

FrontalCortex <- FindVariableFeatures(object = FrontalCortex)

all_cells <- rownames(FrontalCortex)
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 1)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

new.cluster.ids <- c("Ex1", "Oli1", "Ex2", "Ex3", "Ast", "Oli2", "In1", "In2", "In3",
                     "Ex4", "OPC", "Ex5", "Ex6", "Ex7", "In4", "Ex8", "In6", "Mic",
                     "Ex9", "Ex10", "Ex11", "Ex12", "End/Per")
names(new.cluster.ids) <- levels(FrontalCortex)
FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)


plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_findct.png"), device=)

## Doublet removal
## pK indetification
sweep.res.list_pbmc <- paramSweep_v3(FrontalCortex, PCs = 1:20, sct = FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.res.list_pbmc, GT = FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
bcmvn_pbmc$pK <- as.numeric(as.character(bcmvn_pbmc$pK))

## Doublet proportion estimate
annotations <- FrontalCortex@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_pbmc@meta.data$ClusteringResults
nExp_poi <- round(0.15*nrow(FrontalCortex@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

FrontalCortex <- doubletFinder_v3(FrontalCortex, PCs = 1:20, pN = 0.25, pK = bcmvn_pbmc$pK[which.max(bcmvn_pbmc$BCmetric)], nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

DimPlot(FrontalCortex,group.by = colnames(FrontalCortex@meta.data)[grep("DF", colnames(FrontalCortex@meta.data))])
ggsave("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_findct_doublets.png", width = 7, height = 7)

ASD <- SubsetData(ASD, cells=rownames(ASD@meta.data)[which(ASD@meta.data$DF.classification == "Singlet")])
DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
ggsave("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_findct_doubletsremoved.png", width = 8, height = 7)

saveRDS(FrontalCortex, file = path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat_relabeled.rds"))




# FrontalCortex.markers <- FindAllMarkers(FrontalCortex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# FrontalCortex.markers %>% group_by(cluster)
# 
# marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
# all_known_marker_genes = marker_gene_table$Human.Gene
# 
# intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)
# 
# dotplot <- DotPlot(FrontalCortex, features = intersection) + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
# ggsave(path.expand("~/Lake/FrontalCortex/Dotplot_GSE97930_FrontalCortex_Seurat.png"), device=)
# 
# 
# 
# brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
# 
# diff_expressed = FrontalCortex.markers
# 
# celltypelist = brain_genes$Cell.type
# 
# celltypes <- character()
# for (gene in diff_expressed$gene) {
#   if (gene %in% brain_genes$Human.Gene) {
#     celltype = brain_genes$Cell.type[brain_genes$Human.Gene == gene]
#     celltype = paste(celltype, collapse = ', ')
#     celltypes <- c(celltypes, celltype)
#   }
#   else {
#     celltypes <- c(celltypes, "unknown")
#   }
# }
# 
# diff_expressed$cell_type <- celltypes
# 
# diff_expressed_condensed = diff_expressed[diff_expressed$cell_type != "unknown",]
# 
# 
# write.table(diff_expressed, file = path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_diffexpressed.txt"), sep = '\t')
# 
# write.table(diff_expressed_condensed, file = path.expand("~/Lake/FrontalCortex/umap_GSE97930_FrontalCortex_Seurat_diffexpressed_condensed.txt"), sep = '\t')
# 
# FeaturePlot(FrontalCortex, features = c("APBB1IP", "P2RY12"))
# ggsave(path.expand("~/Lake/FrontalCortex/GSE97930_FrontalCortex_Seurat_MicrogliaFeatures.png"), device=, width = 14, height = 7)
