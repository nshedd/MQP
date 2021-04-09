library(dplyr)
library(Seurat)
library(ggplot2)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_SingleR.RDS')

DefaultAssay(BA4.6) <- "RNA"
markers <- FindConservedMarkers(BA4.6, grouping.var = "Group", verbose = FALSE)
write.table(markers, '/data/rusers/sheddn/UCLA-ASD/data/BA4.6_DEGs_byGroup.txt')

markers %>% group_by(cluster) %>% top_n(n = 3, wt = CTL_avg_logFC)
head(markers)



BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_SingleR.RDS')

DefaultAssay(BA9) <- "RNA"
markers <- FindConservedMarkers(BA9, grouping.var = "Group", verbose = FALSE)
write.table(markers, '/data/rusers/sheddn/UCLA-ASD/data/BA9_DEGs_byGroup.txt')

markers %>% group_by(cluster) %>% top_n(n = 3, wt = CTL_avg_logFC)
head(markers)
