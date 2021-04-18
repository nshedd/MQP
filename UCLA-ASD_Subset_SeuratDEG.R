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

clusters_BA4.6 <- c('Ex1','Ex3','In1','OPC','In2','Oli','Ex5','Ex6',
                    'In3','Ast1','In4','Ex7','Ex8','Ex9','In5','Ex10','Ex11','Ast2',
                    'Ex12','Ex13','In6','In7','In8','Ex14','Ex15','Ex16')
                    #No DEGs pass threshold - 'Ex2','Ex4','Ex16'

BA4.6_full <- BA4.6

for (i in clusters_BA4.6) {
  BA4.6_full <- BA4.6
  
  print(i)
  
  sub <- subset(BA4.6_full, idents = i)
  Idents(sub) <- "Group"
  avg.sub <- log1p(AverageExpression(sub, verbose = FALSE)$RNA)
  avg.sub$gene <- rownames(avg.sub)
  
  print(avg.sub[1:10,])
  
  avg.sub$diff = abs(avg.sub$ASD - avg.sub$CTL)
  avg.sub <- avg.sub[order(avg.sub$diff, decreasing=TRUE),] 
  
  genes.to.label <- row.names(avg.sub[1:10,])
  
  p1 <- ggplot(avg.sub, aes(CTL, ASD)) + geom_point() + ggtitle(i)
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/subset/Seurat_integrated/plots/BA4.6_DEGS_",i,'.png', sep='')
  ggsave(link)
}


