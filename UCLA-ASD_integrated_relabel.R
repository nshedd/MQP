library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)


Lake <- readRDS("~/Lake/FrontalCortex/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_Seurat_relabeled.rds")

Lake_SCE <- as.SingleCellExperiment(Lake)
Lake_labels <- Idents(Lake)

# BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')
# 
# print('Running SingleR by cluster...')
# BA4.6_SingleR <- SingleR(test=GetAssayData(BA4.6, assay = "RNA"),
#                          ref=Lake_SCE,
#                          method='cluster',
#                          labels=Lake_labels,
#                          clusters=Idents(BA4.6),
#                          assay.type.test = "logcounts",
#                          assay.type.ref = "logcounts")
# 
# print(BA4.6_SingleR$labels)
# write.table(BA4.6_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/data/BA4.6_cluster_ids.txt")
# 
# new.cluster.ids <- BA4.6_SingleR$labels
# names(new.cluster.ids) <- levels(BA4.6)
# BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)
# 
# print('Plotting by cluster...')
# 
# DimPlot(BA4.6, label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel.png', width = 8, height = 7)
# 
# # new_ids = c("Ex1","Ex2","Oli1","Ast1","Ex3","Ex4","In1","In2","OPC1","Ex5","In3","In4","Ex6","Ex7","Ast2","Ex8","Ex9","Ex10","Ex11","Ex12",
# # "Ex13","Oli2","Ex14","In5","Ast3","In6","Ast4","Ex15","In7","Ex16","Ex17","In8","OPC2","In8","Ex18","OPC3","In9","Ast5","OPC4")
# # 
# # new.cluster.ids <- new_ids
# # names(new.cluster.ids) <- levels(BA4.6)
# # BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)
# 
# DimPlot(BA4.6, reduction = "umap", split.by = "Group")
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel_bycluster_separate.png', width = 16, height = 7)
# 
# print('Running SingleR by cell...')
# BA4.6_SingleR <- SingleR(test=GetAssayData(BA4.6, assay = "RNA"),
#                          ref=Lake_SCE,
#                          labels=Lake_labels,
#                          clusters=Idents(BA4.6),
#                          assay.type.test = "logcounts",
#                          assay.type.ref = "logcounts")
# 
# 
# BA4.6$SingleR.pruned.calls <- BA4.6_SingleR$pruned.labels
# BA4.6$SingleR.calls <- BA4.6_SingleR$labels
# 
# print('Plotting by cell...')
# 
# DimPlot(BA4.6, group.by="SingleR.calls", label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel_bycell.png', width = 8, height = 7)
BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_SingleR.RDS')

DimPlot(BA4.6, group.by="SingleR.calls", reduction = "umap", split.by = "Group", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel_bycell_separate.png', width = 16, height = 7)

# saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_SingleR.RDS')




BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

# print('Running SingleR by cluster...')
# BA9_SingleR <- SingleR(test=GetAssayData(BA9, assay = "RNA"),
#                        ref=Lake_SCE,
#                        method='cluster',
#                        labels=Lake_labels,
#                        clusters=Idents(BA9),
#                        assay.type.test = "logcounts",
#                        assay.type.ref = "logcounts")
# 
# print(BA9_SingleR$labels)
# write.table(BA9_SingleR$labels, "/data/rusers/sheddn/UCLA-ASD/data/BA9_cluster_ids.txt")
# 
# new.cluster.ids <- BA9_SingleR$labels
# names(new.cluster.ids) <- levels(BA9)
# BA9 <- RenameIdents(BA9, new.cluster.ids)
# 
# print('Plotting...')
# 
# DimPlot(BA9, label=TRUE, pt.size=0.5)
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel.png', width = 8, height = 7)
# 
# # new_ids = c("Ex1","Ex2","Oli1","Ast1","Ex3","Ex4","In1","In2","OPC1","Ex5","In3","In4","Ex6","Ex7","Ast2","Ex8","Ex9","Ex10","Ex11","Ex12",
# # "Ex13","Oli2","Ex14","In5","Ast3","In6","Ast4","Ex15","In7","Ex16","Ex17","In8","OPC2","In8","Ex18","OPC3","In9","Ast5","OPC4")
# # 
# # new.cluster.ids <- new_ids
# # names(new.cluster.ids) <- levels(BA9)
# # BA9 <- RenameIdents(BA9, new.cluster.ids)
# 
# DimPlot(BA9, reduction = "umap", split.by = "Group")
# ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel_bycluster_separate.png', width = 16, height = 7)

print('Running SingleR by cell...')
BA9_SingleR <- SingleR(test=GetAssayData(BA9, assay = "RNA"),
                       ref=Lake_SCE,
                       labels=Lake_labels,
                       clusters=Idents(BA9),
                       assay.type.test = "logcounts",
                       assay.type.ref = "logcounts")


BA9$SingleR.pruned.calls <- BA9_SingleR$pruned.labels
BA9$SingleR.calls <- BA9_SingleR$labels

print('Plotting...')

DimPlot(BA9, group.by="SingleR.calls", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel_bycell.png', width = 8, height = 7)

DimPlot(BA9, group.by="SingleR.calls", reduction = "umap", split.by = "Group", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel_bycell_separate.png', width = 16, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_SingleR.RDS')
