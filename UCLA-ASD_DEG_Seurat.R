library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

new_ids = c("Ex1","Ex2","Oli1","Ast1","Ex3","Ex4","In1","In2","OPC1","Ex5","In3","In4","Ex6","Ex7","Ast2","Ex8","Ex9","Ex10","Ex11","Ex12",
            "Ex13","Oli2","Ex14","In5","Ast3","In6","Ast4","Ex15","In7","Ex16","Ex17","In8","OPC2","In8","Ex18","OPC3","In9","Ast5","OPC4")
new.cluster.ids <- new_ids
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

print(Idents(BA4.6))

DefaultAssay(BA4.6) <- "RNA"
markers <- FindConservedMarkers(BA4.6, ident.1="Ex1", grouping.var = "Group", verbose = FALSE)
write.table(markers, '/data/rusers/sheddn/UCLA-ASD/data/BA4.6_DEGs_Ex1_byGroup.txt')

markers %>% group_by(cluster) %>% top_n(n = 3, wt = CTL_avg_logFC)
head(markers)

q()

BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_SingleR.RDS')

DefaultAssay(BA9) <- "RNA"
markers <- FindConservedMarkers(BA9, ident.1="0", grouping.var = "Group", verbose = FALSE)
write.table(markers, '/data/rusers/sheddn/UCLA-ASD/data/BA9_DEGs_byGroup.txt')

markers %>% group_by(cluster) %>% top_n(n = 3, wt = CTL_avg_logFC)
head(markers)
