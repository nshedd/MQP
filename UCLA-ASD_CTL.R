library(dplyr)
library(Seurat)
library(ggplot2)


dir = '/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/11BW_s14/plots'

data_dir <- '/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/Solo.out/Gene/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
CTL = CreateSeuratObject(counts = expression_matrix)


CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")
CTL_VlnPlot <- VlnPlot(CTL, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/qcvlnplot_CTL.png"), device=, width = 14, height = 7)

plot1 <- FeatureScatter(CTL, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(CTL, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot <- plot1 + plot2
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/featurescatter_CTL.png"), device=, width = 14, height = 7)

CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)

CTL <- FindVariableFeatures(object = CTL)

all_cells <- rownames(CTL)
CTL <- ScaleData(CTL, features = all_cells)

CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))

pcs_plot_CTL <- VizDimLoadings(CTL, dims = 1:2, reduction = "pca")
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/qcdimloadings_CTL.png"), device=, width = 14, height = 7)

pca_plot_CTL <- DimPlot(CTL, reduction = "pca")
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/pcascatter_CTL.png"), device=, width = 14, height = 7)

pca_heatmap_CTL <- DimHeatmap(CTL, dims = 1, balanced = TRUE)
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/pcaheatmap_CTL.png"), device=, width = 14, height = 7)

CTL <- JackStraw(CTL, num.replicate = 100)
CTL <- ScoreJackStraw(CTL, dims = 1:20)
Jackstraw_CTL <- JackStrawPlot(CTL, dims = 1:15)
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/jackstrawplot_CTL.png"), device=, width = 14, height = 7)

Elbow_CTL <- ElbowPlot(CTL)
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/elbowplot_CTL.png"), device=, width = 14, height = 7)

CTL <- FindNeighbors(CTL, dims = 1:10)
CTL <- FindClusters(CTL, resolution = 0.5)

CTL <- RunUMAP(CTL, dims = 1:10)

#new.cluster.ids <- c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Astrocyte", "Neuron 5", "? 1", "Oligodendrocyte", "Microglia", "Neuron 6", "? 2", "? 3")
#names(new.cluster.ids) <- levels(CTL)
#CTL <- RenameIdents(CTL, new.cluster.ids)

plot = DimPlot(CTL, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/umap_CTL_orgident.png"), device=)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(CTL.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(CTL, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/Dotplot_CTL.png"), device=)



brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = CTL.markers

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


write.table(diff_expressed, file = path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/umap_CTL_diffexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/umap_CTL_diffexpressed_condensed.txt"), sep = '\t')

FeaturePlot(CTL, features = c("SLC17A7", "GRIN1", "GRIN2B"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_ExcitatoryFeatures.png"), device=)

FeaturePlot(CTL, features = c("GAD1", "GAD2", "SLC6A1"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_InhibitoryFeatures.png"), device=)

FeaturePlot(CTL, features = c("CLDN11", "MOG", "MOBP", "MBP"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_OligodendrocyteFeatures.png"), device=)

FeaturePlot(CTL, features = c("GFAP", "SLC1A2", "SLC1A3", "SLC4A4", "AQP4"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_AstrocyteFeatures.png"), device=)

FeaturePlot(CTL, features = c("COBLL1", "DUSP1", "FLT1"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_EndothelialFeatures.png"), device=)

FeaturePlot(CTL, features = c("COBLL1", "PDGFRB", "COBLL1", "PDGFRB"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_PericyteFeatures.png"), device=)

FeaturePlot(CTL, features = c("PCDH15"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_OPCFeatures.png"), device=)

FeaturePlot(CTL, features = c("APBB1IP", "P2RY12", "PTPRC"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_MicrogliaFeatures.png"), device=)

FeaturePlot(CTL, features = c("RYR1"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_PurkinjeFeatures.png"), device=)

FeaturePlot(CTL, features = c("RELN", "GRM4", "RBFOX3"))
ggsave(path.expand("/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq/10BW_s13/plots/CTL_GranuleFeatures.png"), device=)
