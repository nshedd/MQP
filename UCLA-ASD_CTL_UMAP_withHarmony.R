library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)

## BA4/6
# CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6')
# 
# CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")
# 
# CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(CTL)
# CTL <- ScaleData(CTL, features = all.genes)

# CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))
# 
# CTL <- RunHarmony(CTL, "orig.ident")
# 
# print("Saving Harmony data BA4.6...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA4.6.RDS')

print("Loading Harmony data BA4.6...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA4.6.RDS')

CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)

CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution=0.5)

DimPlot(CTL, group.by="ident")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_lowerres.png', width = 8, height = 7)

print("Saving UMAP data BA4.6...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_lowerres.RDS')

# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6.RDS')

# new.cluster.ids <- c('Ex1','Ex2','Ex3','Ex4','Oli','Ex5','In1','In2','In3','OPC','Ex6',
#                     'Ex7','Ast1','End?','Ex8','Mic','Ex9','In4','Ex10','In5','Ex11','In6','Per?',
#                      'In6','Dop?','In7','Ex12','Ast2','Ex13','Ast3','Ex14','Ex15','Ex16','End/Per')
# names(new.cluster.ids) <- levels(CTL)
# CTL <- RenameIdents(CTL, new.cluster.ids)
# 
# DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA4.6_labeled.png', width = 8, height = 7)

CTL.markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CTL.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(CTL.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(CTL, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/CTL_BA4.6_dotplot_lowerres.png", width = 14, height = 7)


## BA9
# CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_BA9')
# 
# CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")
# 
# CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)
# 
# CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)
# 
# all.genes <- rownames(CTL)
# CTL <- ScaleData(CTL, features = all.genes)
# 
# CTL <- RunPCA(CTL, features = VariableFeatures(object = CTL))
# 
# CTL <- RunHarmony(CTL, "orig.ident", max.iter.harmony=20)
# 
# print("Saving Harmony data BA9...")
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA9.RDS')

print("Loading Harmony data BA9...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BA9.RDS')

CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)

CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution=0.5)

DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_BA9_lowerres.png', width = 8, height = 7)

print("Saving UMAP data BA9...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_lowerres.RDS')

CTL.markers <- FindAllMarkers(CTL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CTL.markers %>% group_by(cluster)

marker_gene_table = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")
all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(CTL.markers$gene, all_known_marker_genes)

dotplot <- DotPlot(CTL, features = intersection) + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_y_discrete(limits = rev(levels(CTL$seurat_clusters)))
ggsave("/data/rusers/sheddn/UCLA-ASD/plots/CTL_BA9_dotplot_lowerres.png", width = 14, height = 7)



# CTL = readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_SampleLabels')

# CTL[["percent.mt"]] <- PercentageFeatureSet(CTL, pattern = "^MT-")

# CTL <- NormalizeData(CTL, normalization.method = "LogNormalize", scale.factor = 10000)

# CTL <- FindVariableFeatures(CTL, selection.method = "vst", nfeatures = 2000)

# all.genes <- rownames(CTL)
# CTL <- ScaleData(CTL, features = all.genes)

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_ScaledData.RDS')

# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_PCAprocessed_BySample.RDS')

# CTL <- RunHarmony(CTL, "orig.ident")

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BySample.RDS')

# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BySample.RDS')

# CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)

# CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters()

# DimPlot(CTL, group.by="ident")
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony_samples.png', width = 8, height = 7)

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_OrigIdent.RDS')

# CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_Harmonyprocessed_BySample.RDS')

# CTL <- RunUMAP(CTL, reduction = "harmony", dims = 1:30)

# CTL <- FindNeighbors(CTL, reduction = "harmony", dims = 1:30) %>% FindClusters()

# DimPlot(CTL, group.by="ident", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CTL-UMAP_Harmony.png', width = 8, height = 7)

# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony.RDS')
