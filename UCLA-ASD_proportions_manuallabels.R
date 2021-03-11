library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)

## ASD-BA4.6
print("Loading ASD UMAP data w/o Doublets BA4.6...")
ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved_Relabeled.RDS')

new.cluster.ids <- c('Ex','Ex','Ex','Ast','Oli','In','Ex','In','Ex','Ex','In',
                     'Ex','OPC','End/Per','Ex','Mic','In','Ex','Ex','Ex','Ex','Ex',
                     'Oli','Ex','Ex','Ast','In','In','Ast','Mic','In','Ex','In',
                     'In','In','Ast','Mic','In','Ex','Ast','Ast','Ast','OPC','Mic')
names(new.cluster.ids) <- levels(ASD)
ASD <- RenameIdents(ASD, new.cluster.ids)

ASD_BA4.6_prop_manual_list = Idents(ASD)
ASD_BA4.6_prop_manual_table = table(ASD_BA4.6_prop_manual_list)

ASD_BA4.6_prop_df = ASD_BA4.6_prop_manual_table %>% as.data.frame
colnames(ASD_BA4.6_prop_df) <- c('Cell_Type', 'FreqASD')
print(ASD_BA4.6_prop_df)


## CTL-BA4.6
print("Loading CTL UMAP data w/o Doublets BA4.6...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved_Relabeled.RDS')

new.cluster.ids <- c('Ex','In','Ex','Ex','Ex','Ex','Oli','In','End/Per',
                    'Ex','In','Ex','OPC','Ast','In','Ex','Mic','Ex','Ex',
                    'Ex','In','In','Ex','Ex','In','Ex','Ast','Ex')
names(new.cluster.ids) <- levels(CTL)
CTL <- RenameIdents(CTL, new.cluster.ids)

CTL_BA4.6_prop_manual_list = Idents(CTL)
CTL_BA4.6_prop_manual_table = table(CTL_BA4.6_prop_manual_list)

CTL_BA4.6_prop_df = CTL_BA4.6_prop_manual_table %>% as.data.frame
colnames(CTL_BA4.6_prop_df) <- c('Cell_Type', 'FreqCTL')
print(CTL_BA4.6_prop_df)


BA4.6_prop_df = merge(x = ASD_BA4.6_prop_df, y = CTL_BA4.6_prop_df, by = 'Cell_Type', all = TRUE)
print(BA4.6_prop_df)

BA4.6_prop_df = melt(BA4.6_prop_df, id.vars='Cell_Type')

ggplot(data=BA4.6_prop_df, aes(x = Cell_Type, y = value, fill=variable)) + geom_bar(stat="identity") +ggtitle("BA 4/6 - Cluster Labels")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CellTypeProportion_BA4.6.png')


## ASD-BA9
## cuz I'm a dummy and forgot to save these files 
print("Loading UMAP data w/o Doublets BA9...")
ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved.RDS')

new.cluster.ids <- c('Ex1','Ex2','Ex3','In1','OPC','Ex4','Ast1','In2','Ex5',
                     'Mic','Ex6','Ex7','In3','Ex8','Ex9','Ex10','Oli1','Ast2','In4',
                     'Ex11','End/Per1','Ex12','In5','Ex13','Ex14','End/Per2','Oli2','Oli3')
names(new.cluster.ids) <- levels(ASD)
ASD <- RenameIdents(ASD, new.cluster.ids)

print("Saving UMAP data w/o Doublets BA9...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved_Relabeled.RDS')

new.cluster.ids <- c('Ex','Ex','Ex','In','OPC','Ex','Ast','In','Ex',
                     'Mic','Ex','Ex','In','Ex','Ex','Ex','Oli','Ast','In',
                     'Ex','End/Per','Ex','In','Ex','Ex','End/Per','Oli','Oli')
names(new.cluster.ids) <- levels(ASD)
ASD <- RenameIdents(ASD, new.cluster.ids)

DimPlot(ASD, group.by="ident", label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/ASD-UMAP_Harmony_BA9_labeled.png', width = 8, height = 7)

ASD_BA9_prop_manual_list = Idents(ASD)
ASD_BA9_prop_manual_table = table(ASD_BA9_prop_manual_list)

ASD_BA9_prop_df = ASD_BA9_prop_manual_table %>% as.data.frame
colnames(ASD_BA9_prop_df) <- c('Cell_Type', 'FreqASD')


## CTL-BA9
print("Saving UMAP data w/o Doublets BA9...")
CTL <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved.RDS')

new.cluster.ids <- c('Ex1','Ex2','Ex3','Ast1','OPC','In1','Ex4','Ex5','In2', 
                     'Ex6','In3','Ex7','Ast2','Mic','In4','End/Per1','Ex8','Oli1',
                     'Ex9','Ex10','In5','In6','End/Per2','Ex11','In7','Ex12','Oli2')
names(new.cluster.ids) <- levels(CTL)
CTL <- RenameIdents(CTL, new.cluster.ids)

print("Saving UMAP data w/o Doublets BA9...")
saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA9_DoubletsRemoved_Relabeled.RDS')

new.cluster.ids <- c('Ex','Ex','Ex','Ast','OPC','In','Ex','Ex','In', 
                     'Ex','In','Ex','Ast','Mic','In','End/Per','Ex','Oli',
                     'Ex','Ex','In','In','End/Per','Ex','In','Ex','Oli')
names(new.cluster.ids) <- levels(CTL)
CTL <- RenameIdents(CTL, new.cluster.ids)

CTL_BA9_prop_manual_list = Idents(CTL)
CTL_BA9_prop_manual_table = table(CTL_BA9_prop_manual_list)

CTL_BA9_prop_df = CTL_BA9_prop_manual_table %>% as.data.frame
colnames(CTL_BA9_prop_df) <- c('Cell_Type', 'FreqCTL')


BA9_prop_df = melt(BA9_prop_df, id.vars='Cell_Type')

ggplot(data=BA9_prop_df, aes(x = Cell_Type, y = value, fill=variable)) + geom_bar(stat="identity") +ggtitle("BA 9 - Cluster Labels")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/CellTypeProportion_BA9.png')
