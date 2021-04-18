library(Seurat)
library(ggplot2)
library(cowplot)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/data/BA4.6_WithDEGs.RDS')

clusters_BA4.6 <- c('Ex1','Ex2','Ex3','In1','Ex4','OPC','In2','Oli','Ex5','Ex6',
                    'In3','Ast1','In4','Ex7','Ex8','Ex9','In5','Ex10','Ex11','Ast2',
                    'Ex12','Ex13','In6','In7','In8','Ex14','Ex15','Ex16')

new.cluster.ids <- clusters_BA4.6
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

for (i in clusters_BA4.6) {
  sub <- subset(BA4.6, idents = i)
  Idents(sub) <- "Group"
  print(Idents(sub))
  avg.sub <- log1p(AverageExpression(sub, verbose = FALSE)$RNA)
  avg.sub$gene <- rownames(avg.sub)
  
  Idents(BA4.6) <- "ident.Group"
  print(Idents(BA4.6))
  interferon.response <- FindMarkers(BA4.6, ident.1 = paste(i, "ASD", sep="."), ident.2 = paste(i, "CTL", sep="."), verbose = FALSE)
  
  genes.to.label <- row.names(interferon.response[1:10,])
  
  p1 <- ggplot(avg.sub, aes(ASD, CTL)) + geom_point() + ggtitle(i)
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/BA4.6_DEGS_",i,'.png', sep='')
  ggsave(link)
}


