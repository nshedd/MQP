library(dplyr)
library(Seurat)
library(ggplot2)


dir = '/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots'

data_dir <- '/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/Solo.out/Gene/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
ASD = CreateSeuratObject(counts = expression_matrix)


ASD[["percent.mt"]] <- PercentageFeatureSet(ASD, pattern = "^MT-")
ASD_VlnPlot <- VlnPlot(ASD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path.expand(dir"/qcvlnplot_ASD.png"), device=, width = 14, height = 7)

plot1 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ASD, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave(path.expand(dir"/featurescatter_ASD.png"), device=, width = 14, height = 7)

ASD <- NormalizeData(ASD, normalization.method = "LogNormalize", scale.factor = 10000)

ASD <- FindVariableFeatures(object = ASD)

all_cells <- rownames(ASD)
ASD <- ScaleData(ASD, features = all_cells)

ASD <- RunPCA(ASD, features = VariableFeatures(object = ASD))

pcs_plot_ASD <- VizDimLoadings(ASD, dims = 1:2, reduction = "pca")
ggsave(path.expand(dir"/qcdimloadings_ASD.png"), device=, width = 14, height = 7)

pca_plot_ASD <- DimPlot(ASD, reduction = "pca")
ggsave(path.expand(dir"/pcascatter_ASD.png"), device=, width = 14, height = 7)

pca_heatmap_ASD <- DimHeatmap(ASD, dims = 1, balanced = TRUE)
ggsave(path.expand(dir"/pcaheatmap_ASD.png"), device=, width = 14, height = 7)

ASD <- JackStraw(ASD, num.replicate = 100)
ASD <- ScoreJackStraw(ASD, dims = 1:20)
Jackstraw_ASD <- JackStrawPlot(ASD, dims = 1:15)
ggsave(path.expand(dir"/jackstrawplot_ASD.png"), device=, width = 14, height = 7)

Elbow_ASD <- ElbowPlot(ASD)
ggsave(path.expand(dir"/elbowplot_ASD.png"), device=, width = 14, height = 7)

ASD <- FindNeighbors(ASD, dims = 1:10)
ASD <- FindClusters(ASD, resolution = 0.5)

ASD <- RunUMAP(ASD, dims = 1:10)

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(ASD)
#ASD <- RenameIdents(ASD, new.cluster.ids)

plot = DimPlot(ASD, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand(dir"/umap_ASD_orgident.png"), device=)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(ASD.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(ASD, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(ASD$seurat_clusters)))
ggsave(path.expand(dir"/Dotplot_ASD.png"), device=)



brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = ASD.markers

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


write.table(diff_expressed, file = path.expand(dir"/umap_ASD_diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand(dir"/umap_ASD_diffexpressed_condensed.txt"), sep = '\t')

FeaturePlot(ASD, features = c("SLC17A7", "GRIN1", "GRIN2B"))
ggsave(path.expand(dir"/ASD_ExcitatoryFeatures.png"), device=)

FeaturePlot(ASD, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand(dir"/ASD_InhibitoryFeatures.png"), device=)

FeaturePlot(ASD, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand(dir"/ASD_OligodendrocyteFeatures.png"), device=)

FeaturePlot(ASD, features = c("GFAP", "SLC1A2", "SLC1A3", "SLC4A4", "AQP4"))
ggsave(path.expand(dir"/ASD_AstrocyteFeatures.png"), device=)

FeaturePlot(ASD, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand(dir"/ASD_EndothelialFeatures.png"), device=)

FeaturePlot(ASD, features = c("COBLL1", "PDGFRB", "COBLL1", "PDGFRB"))
ggsave(path.expand(dir"/ASD_PericyteFeatures.png"), device=)

FeaturePlot(ASD, features = c("PCDH15"))
ggsave(path.expand(dir"/ASD_OPCFeatures.png"), device=)

FeaturePlot(ASD, features = c("APBB1IP", "P2RY12", "PTPRC"))
ggsave(path.expand(dir"/ASD_MicrogliaFeatures.png"), device=)

FeaturePlot(ASD, features = c("RYR1"))
ggsave(path.expand(dir"/ASD_PurkinjeFeatures.png"), device=)

FeaturePlot(ASD, features = c("RELN", "GRM4", "RBFOX3"))
ggsave(path.expand(dir"/ASD_GranuleFeatures.png"), device=)
