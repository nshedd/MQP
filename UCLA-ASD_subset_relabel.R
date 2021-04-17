library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)


BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6_WithDEGs.RDS')

clusters_BA4.6 <- c('Ex1','Ex2','Ex3','In1','In2','Ast1','In3','Ex4','Ex5','In4',
                    'Oli1','Ex6','Ex7','Ex8','Ex9','Ex10','OPC1','Ast2','In5','Ex11',
                    'Ex12','Ex13','Ex14','In6','In7','Ast3','OPC2','Ex15','In8','Oli2',
                    'In8','Ex16','Ex17','In9','Ast4','In10','Ast5','Ex18','Ex19','Ast6',
                    'In11','Ex20','Ex21','OPC3','In12','Ex22','Ex23','Ex24','In13','OPC4','In14','Ex25')

new.cluster.ids <- clusters_BA4.6
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

DimPlot(BA4.6, label=TRUE, pt.size=0.5, split.by="Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel.png', width = 16, height = 7)

saveRDS(BA4.6, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6_Relabeled.RDS')


BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9_WithDEGs.RDS')

clusters_BA9 <- c('Ex1','Ex2','Ex3','In1','In2','Ex4','Ex5','Ast1','In3','In4',
                  'In5','In6','Ex6','Ex7','Ex8','Ex9','Ex10','Ex11','Ex12','Ex13',
                  'Ex14','Ex15','Ast2','OPC1','Ex16','Ex17','Ex18','In7','In8','Ex19',
                  'Ex20','Ast3','Ex21','In9','In10','OPC2','In11','Ex22','In12','Ex23')

new.cluster.ids <- clusters_BA9
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

DimPlot(BA9, label=TRUE, pt.size=0.5, split.by="Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/subset/plots/UMAP_Harmony_BA9_integrated_SingleRlabel.png', width = 16, height = 7)

saveRDS(BA9, '/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9_Relabeled.RDS')
