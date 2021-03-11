library(dplyr)
library(Seurat)
library(ggplot2)
library(harmony)
library(DoubletFinder)
library(SingleR)
library(SingleCellExperiment)

## ASD-BA4.6
## cuz I'm a dummy and forgot to save these files 
print("Loading UMAP data w/o Doublets BA4.6...")
ASD <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved.RDS')

new.cluster.ids <- c('Ex1','Ex2','Ex3','Ast1','Oli1','In1','Ex4','In2','Ex5','Ex6','In3',
                     'Ex7','OPC1','End/Per','Ex8','Mic1','In4','Ex9','Ex10','Ex11','Ex12','Ex13',
                     'Oli2','Ex14','Ex15','Ast2','In5','In6','Ast3','Mic2','In7','Ex16','In8',
                     'In9','In10','Ast4','Mic3','In11','Ex17','Ast5','Ast6','Ast7','OPC2','Mic4')
names(new.cluster.ids) <- levels(ASD)
ASD <- RenameIdents(ASD, new.cluster.ids)

print("Saving UMAP data w/o Doublets BA4.6...")
saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/CTL_UMAPprocessed_BySample_Harmony_BA4.6_DoubletsRemoved_Relabeled.RDS')

new.cluster.ids <- c('Ex','Ex','Ex','Ast','Oli','In','Ex','In','Ex','Ex','In',
                     'Ex','OPC','End/Per','Ex','Mic','In','Ex','Ex','Ex','Ex','Ex',
                     'Oli','Ex','Ex','Ast','In','In','Ast','Mic','In','Ex','In',
                     'In','In','Ast','Mic','In','Ex','Ast','Ast','Ast','OPC','Mic')
names(new.cluster.ids) <- levels(ASD)
ASD <- RenameIdents(ASD, new.cluster.ids)

ASD_BA4.6_prop_manual = Idents(ASD)
ASD_BA4.6_prop_manual = table(ASD_BA4.6_prop)
print(ASD_BA4.6_prop)


