library(ggplot2)
library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_DoubletsRemoved')

clusters_BA4.6 = c("Ex1","Ex2","Oli1","Ast1","Ex3","Ex4","In1","In2","OPC1",
            "Ex5","In3","In4","Ex6","Ex7","Ast2","Ex8","Ex9","Ex10","In5",
            "Ex12","Ex13","Ast2","Ex14","In6","Ast3","In7","Ast4","Mic","In8",
            "Ex15","Ex16","In9","OPC2","In10","Ex17","OPC3","In11","Ast5","OPC4")

new.cluster.ids <- clusters_BA4.6
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

clusters_BA4.6 = c("Ex1","Ex2","Oli1","Ast1","Ex3","Ex4","In1","In2","OPC1",
            "Ex5","In3","In4","Ex6","Ex7","Ast2","Ex8","Ex9","Ex10","In5",
            "Ex12","Ex13","Ast2","Ex14","In6","Ast3","In7","Ast4","Mic","In8",
            "Ex15","Ex16","In9","OPC2","In10","Ex17","OPC3","In11","Ast5","OPC4")
            #No DEGs pass threshold -

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
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG/BA4.6_DEGS_",i,'.png', sep='')
  ggsave(link)
}


BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_DoubletsRemoved')

clusters_BA9 = c('Ex1','Ex2','Ex3','In1','Ast1','OPC1','Ex4','Mic','In2','Ex5',
            'Ex6','In3','Ex7','Ex8','Ex9','Ast2','Oli1','Ex10','In4','Ex11',
            'Ex12','Ex13','In5','Ex14','Ex15','Ex16','In6','In7','Ex17','In8','Ex18')

new.cluster.ids <- clusters_BA9
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

clusters_BA9 = c('Ex1','Ex2','Ex3','In1','Ast1','OPC1','Ex4','Mic','In2','Ex5',
            'Ex6','In3','Ex7','Ex8','Ex9','Ast2','Oli1','Ex10','In4','Ex11',
            'Ex12','Ex13','In5','Ex14','Ex15','Ex16','In6','In7','Ex17','In8','Ex18')
            #No DEGs pass threshold -

BA9_full <- BA9

for (i in clusters_BA9) {
  BA9_full <- BA9
  
  print(i)
  
  sub <- subset(BA9_full, idents = i)
  Idents(sub) <- "Group"
  avg.sub <- log1p(AverageExpression(sub, verbose = FALSE)$RNA)
  avg.sub$gene <- rownames(avg.sub)
  
  print(avg.sub[1:10,])
  
  avg.sub$diff = abs(avg.sub$ASD - avg.sub$CTL)
  avg.sub <- avg.sub[order(avg.sub$diff, decreasing=TRUE),] 
  
  genes.to.label <- row.names(avg.sub[1:10,])
  
  p1 <- ggplot(avg.sub, aes(CTL, ASD)) + geom_point() + ggtitle(i)
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG/BA9_DEGS_",i,'.png', sep='')
  ggsave(link)
}
