library(dplyr)
library(Seurat)
library(ggplot2)


dir = '/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots'

data_dir <- '/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/Solo.out/Gene/filtered'
list.files(data_dir) # Should show barcodes.tsv.gz, features.tsv.gz, and matrix.mtx.gz
data <- Read10X(data.dir = data_dir)
ASD = CreateSeuratObject(counts = data$`Gene Expression`)


FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")
FrontalCortex_VlnPlot <- VlnPlot(FrontalCortex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path.expand("dir/qcvlnplot_ASD.png"), device=, width = 14, height = 7)

plot1 <- FeatureScatter(FrontalCortex, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(FrontalCortex, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave(path.expand("dir/featurescatter_ASD.png"), device=, width = 14, height = 7)

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)

FrontalCortex <- FindVariableFeatures(object = FrontalCortex)

all_cells <- rownames(FrontalCortex)
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))

pcs_plot_FrontalCortex <- VizDimLoadings(FrontalCortex, dims = 1:2, reduction = "pca")
ggsave(path.expand("dir/qcdimloadings_ASD.png"), device=, width = 14, height = 7)

pca_plot_FrontalCortex <- DimPlot(FrontalCortex, reduction = "pca")
ggsave(path.expand("dir/pcascatter_ASD.png"), device=, width = 14, height = 7)

pca_heatmap_FrontalCortex <- DimHeatmap(FrontalCortex, dims = 1, balanced = TRUE)
ggsave(path.expand("dir/pcaheatmap_ASD.png"), device=, width = 14, height = 7)

FrontalCortex <- JackStraw(FrontalCortex, num.replicate = 100)
FrontalCortex <- ScoreJackStraw(FrontalCortex, dims = 1:20)
Jackstraw_FrontalCortex <- JackStrawPlot(FrontalCortex, dims = 1:15)
ggsave(path.expand("dir/jackstrawplot_ASD.png"), device=, width = 14, height = 7)

Elbow_FrontalCortex <- ElbowPlot(FrontalCortex)
ggsave(path.expand("dir/elbowplot_ASD.png"), device=, width = 14, height = 7)

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:10)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 0.5)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:10)

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(FrontalCortex)
#FrontalCortex <- RenameIdents(FrontalCortex, new.cluster.ids)

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("dir/umap_ASD_orgident.png"), device=)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(FrontalCortex.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(FrontalCortex, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(FrontalCortex$seurat_clusters)))
ggsave(path.expand("dir/Dotplot_ASD.png"), device=)



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


write.table(diff_expressed, file = path.expand("dir/umap_ASD_diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("dir/umap_ASD_diffexpressed_condensed.txt"), sep = '\t')

FeaturePlot(CerebellarHem, features = c("SLC17A7", "GRIN1", "GRIN2B"))
ggsave(path.expand("dir/ASD_ExcitatoryFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand("dir/ASD_InhibitoryFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand("dir/ASD_OligodendrocyteFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("GFAP", "SLC1A2", "SLC1A3", "SLC4A4", "AQP4"))
ggsave(path.expand("dir/ASD_AstrocyteFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand("dir/ASD_EndothelialFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("COBLL1", "PDGFRB", "COBLL1", "PDGFRB"))
ggsave(path.expand("dir/ASD_PericyteFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("PCDH15"))
ggsave(path.expand("dir/ASD_OPCFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("APBB1IP", "P2RY12", "PTPRC"))
ggsave(path.expand("dir/ASD_MicrogliaFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("RYR1"))
ggsave(path.expand("dir/ASD_PurkinjeFeatures.png"), device=)

FeaturePlot(CerebellarHem, features = c("RELN", "GRM4", "RBFOX3"))
ggsave(path.expand("dir/ASD_GranuleFeatures.png"), device=)
