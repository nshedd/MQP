library(ggplot2)
library(Seurat)

BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

clusters_BA9 = c("Ex","Ex","Ex","In","Ast","OPC","Mic","Ex","In","Ex",
                 "In","Ex","Ex","In","Ex","Oli","Ast","Ex","Ex","Ex",
                 "Ex","Ex","Ex","Ex","Ast","In","In","Ex","Ex","In",
                 "OPC","Ex")

new.cluster.ids <- clusters_BA9
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

clusters_BA9 = c('Ast','Ex','In','Mic','OPC','Oli')

BA9_full <- BA9

for (i in clusters_BA9) {
  BA9_full <- BA9
  
  print(i)
  
  sub <- subset(BA9_full, idents = i)
  Idents(sub) <- "Group"
  avg.sub <- log1p(AverageExpression(sub, verbose = FALSE)$RNA)
  avg.sub$gene <- rownames(avg.sub)
  
  BA9_full$celltype.group <- paste(Idents(BA9_full),BA9_full$Group, sep = "_")
  BA9_full$celltype <- Idents(BA9_full)
  Idents(BA9_full) <- "celltype.group"
  response <- FindMarkers(BA9_full, ident.1 = paste(i,"ASD", sep='_'), ident.2 = paste(i,"CTL", sep='_'), test.use="DESeq2", verbose = FALSE)
  print(response[1:25,])
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA9_DEGs-FDR_",i,'.txt', sep='')
  write.table(response, link)
  
  genes.to.label <- row.names(response[1:25,])
  
  p1 <- ggplot(avg.sub, aes(CTL, ASD)) + geom_point() + ggtitle(i)
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA9_DEGS_",i,'.png', sep='')
  ggsave(link)
}

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

clusters_BA4.6 = c("Ex","Ex","Ex","In","Ast","Oli","In","OPC","In","Ex",
                   "Ex","Ex","Ex","Ex","Ex","In","Ex","Ex","In","Ex",
                   "Ex","Oli","Ex","Ast","Ast","In","Ast","In","Mic","In",
                   "Oli","In","Ex","In","In","OPC","In","Ex","Ast","In",
                   "Ast","OPC","Mic")

new.cluster.ids <- clusters_BA4.6
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

clusters_BA4.6 = c('Ast','Ex','In','Mic','OPC','Oli')

BA4.6_full <- BA4.6

for (i in clusters_BA4.6) {
  BA4.6_full <- BA4.6
  
  print(i)
  
  sub <- subset(BA4.6_full, idents = i)
  Idents(sub) <- "Group"
  avg.sub <- log1p(AverageExpression(sub, verbose = FALSE)$RNA)
  avg.sub$gene <- rownames(avg.sub)
  
  BA4.6_full$celltype.group <- paste(Idents(BA4.6_full),BA4.6_full$Group, sep = "_")
  BA4.6_full$celltype <- Idents(BA4.6_full)
  Idents(BA4.6_full) <- "celltype.group"
  response <- FindMarkers(BA4.6_full, ident.1 = paste(i,"ASD", sep='_'), ident.2 = paste(i,"CTL", sep='_'), test.use="DESeq2", verbose = FALSE)
  print(response[1:25,])
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_",i,'.txt', sep='')
  write.table(response, link)
  
  genes.to.label <- row.names(response[1:25,])
  
  p1 <- ggplot(avg.sub, aes(CTL, ASD)) + geom_point() + ggtitle(i)
  p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGS_",i,'.png', sep='')
  ggsave(link)
}


q()


clusters_BA4.6 = c("Ex1","Ex2","Ex3","In1","Ast1","Oli1","In2","OPC1","In3","Ex4",
                   "Ex5","Ex6","Ex7","Ex8","Ex9","In4","Ex10","Ex11","In5","Ex12",
                   "Ex13","Oli2","Ex14","Ast2","Ast3","In6","Ast4","In7","Mic1","In8",
                   "Oli3","In9","Ex15","In10","In11","OPC2","In12","Ex16","Ast5","In13",
                   "Ast6","OPC3","Mic2")


new.cluster.ids <- clusters_BA4.6
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

clusters_BA4.6 = c("Ex1","Ex2","Ex3","In1","Ast1","Oli1","In2","OPC1","In3","Ex4",
                   "Ex5","Ex6","Ex7","Ex8","Ex9","In4","Ex10","Ex11","In5","Ex12",
                   "Ex13","Oli2","Ex14","Ast2","Ast3","In6","Ast4","In7","Mic1","In8",
                   "Oli3","In9","Ex15","In10","In11","OPC2","In12","Ex16","Ast5","In13")
# only in 1 group - "Ast6", "OPC3", "Mic2"

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
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGS_",i,'.png', sep='')
  ggsave(link)
}

clusters_BA9 = c("Ex1","Ex2","Ex3","In1","Ast1","OPC1","Mic","Ex4","In2","Ex5",
                 "In3","Ex6","Ex7","In4","Ex8","Oli","Ast2","Ex9","Ex10","Ex11",
                 "Ex11","Ex12","Ex13","Ex14","Ast3","In5","In6","Ex15","Ex16","In7",
                 "OPC2","Ex17")
new.cluster.ids <- clusters_BA9
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

clusters_BA9 = c("Ex1","Ex2","Ex3","In1","Ast1","OPC1","Mic","Ex4","In2","Ex5",
                 "In3","Ex6","Ex7","In4","Ex8","Oli","Ast2","Ex9","Ex10","Ex11",
                 "Ex11","Ex12","Ex13","Ex14","Ast3","In5","In6","Ex15","Ex16","In7",
                 "OPC2","Ex17")

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
  
  link = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA9_DEGS_",i,'.png', sep='')
  ggsave(link)
}

