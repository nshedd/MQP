library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

group = c('CTL','ASD')
sample = unique(BA4.6$Sample)
cluster = unique(Idents(BA4.6))

for (c in cluster) {
  data <- subset(BA4.6, ident=c)
  expression.data <- as.matrix(GetAssayData(data, slot = "counts"))
  expression_file <- paste("/data/rusers/sheddn/UCLA-ASD/plots/BA4.6_averageexpression_cluster", i, ".txt", sep='')
  write.table(expression.data, expression_file, sep='\t')
  
  head(expression.data)
            
  df = data[,"Group", drop=FALSE]
  df$Sample = data$Sample
  meta_file <- paste("/data/rusers/sheddn/UCLA-ASD/plots/BA4.6_metadata_cluster", i, ".txt", sep='')
  write.table(df, meta_file, sep='\t')
  
  head(df)
}


q()

BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

data.expression = AverageExpression(BA9)
data.expression.rna = data.expression[["RNA"]]
head(data.expression.rna)

cluster.expression.0 = data.expression.rna[,"0", drop=FALSE]
cluster.expression.1 = data.expression.rna[,"1", drop=FALSE]
cluster.expression.2 = data.expression.rna[,"2", drop=FALSE]
cluster.expression.3 = data.expression.rna[,"3", drop=FALSE]
cluster.expression.4 = data.expression.rna[,"4", drop=FALSE]
cluster.expression.5 = data.expression.rna[,"5", drop=FALSE]
cluster.expression.6 = data.expression.rna[,"6", drop=FALSE]
cluster.expression.7 = data.expression.rna[,"7", drop=FALSE]
cluster.expression.8 = data.expression.rna[,"8", drop=FALSE]
cluster.expression.9 = data.expression.rna[,"9", drop=FALSE]
cluster.expression.10 = data.expression.rna[,"10", drop=FALSE]
cluster.expression.11 = data.expression.rna[,"11", drop=FALSE]
cluster.expression.12 = data.expression.rna[,"12", drop=FALSE]
cluster.expression.13 = data.expression.rna[,"13", drop=FALSE]
cluster.expression.14 = data.expression.rna[,"14", drop=FALSE]
cluster.expression.15 = data.expression.rna[,"15", drop=FALSE]
cluster.expression.16 = data.expression.rna[,"16", drop=FALSE]
cluster.expression.17 = data.expression.rna[,"17", drop=FALSE]
cluster.expression.18 = data.expression.rna[,"18", drop=FALSE]
cluster.expression.19 = data.expression.rna[,"19", drop=FALSE]
cluster.expression.20 = data.expression.rna[,"20", drop=FALSE]
cluster.expression.21 = data.expression.rna[,"21", drop=FALSE]
cluster.expression.22 = data.expression.rna[,"22", drop=FALSE]
cluster.expression.23 = data.expression.rna[,"23", drop=FALSE]
cluster.expression.24 = data.expression.rna[,"24", drop=FALSE]
cluster.expression.25 = data.expression.rna[,"25", drop=FALSE]
cluster.expression.26 = data.expression.rna[,"26", drop=FALSE]
cluster.expression.27 = data.expression.rna[,"27", drop=FALSE]
cluster.expression.28 = data.expression.rna[,"28", drop=FALSE]
cluster.expression.29 = data.expression.rna[,"29", drop=FALSE]
cluster.expression.30 = data.expression.rna[,"30", drop=FALSE]

write.table(data.expression.rna, "/data/rusers/sheddn/UCLA-ASD/data/BA9_averageexpression.txt")



group = c('CTL','ASD')
sample = unique(BA9$Sample)
cluster = unique(Idents(BA9))
               
for (s in sample) {
  data <- subset(BA9, subset= Sample==s)
  data.expression = AverageExpression(data)
  data.expression.rna = data.expression[["RNA"]]
  cluster = unique(Idents(data))
  n = paste(data$Group[1], s, sep='_')
  for (c in cluster) {
    if (c == "0"){
      cluster.expression.0[,n] = data.expression.rna[,"0"]
    }
    if (c == "1"){
      cluster.expression.1[,n] = data.expression.rna[,"1"]
    }
    if (c == "2"){
      cluster.expression.2[,n] = data.expression.rna[,"2"]
    }
    if (c == "3"){
      cluster.expression.3[,n] = data.expression.rna[,"3"]
    }
    if (c == "4"){
      cluster.expression.4[,n] = data.expression.rna[,"4"]
    }
    if (c == "5"){
      cluster.expression.5[,n] = data.expression.rna[,"5"]
    }
    if (c == "6"){
      cluster.expression.6[,n] = data.expression.rna[,"6"]
    }
    if (c == "7"){
      cluster.expression.7[,n] = data.expression.rna[,"7"]
    }
    if (c == "8"){
      cluster.expression.8[,n] = data.expression.rna[,"8"]
    }
    if (c == "9"){
      cluster.expression.9[,n] = data.expression.rna[,"9"]
    }
    if (c == "10"){
      cluster.expression.10[,n] = data.expression.rna[,"10"]
    }
    if (c == "11"){
      cluster.expression.11[,n] = data.expression.rna[,"11"]
    }
    if (c == "12"){
      cluster.expression.12[,n] = data.expression.rna[,"12"]
    }
    if (c == "13"){
      cluster.expression.13[,n] = data.expression.rna[,"13"]
    }
    if (c == "14"){
      cluster.expression.14[,n] = data.expression.rna[,"14"]
    }
    if (c == "15"){
      cluster.expression.15[,n] = data.expression.rna[,"15"]
    }
    if (c == "16"){
      cluster.expression.16[,n] = data.expression.rna[,"16"]
    }
    if (c == "17"){
      cluster.expression.17[,n] = data.expression.rna[,"17"]
    }
    if (c == "18"){
      cluster.expression.18[,n] = data.expression.rna[,"18"]
    }
    if (c == "19"){
      cluster.expression.19[,n] = data.expression.rna[,"19"]
    }
    if (c == "20"){
      cluster.expression.20[,n] = data.expression.rna[,"20"]
    }
    if (c == "21"){
      cluster.expression.21[,n] = data.expression.rna[,"21"]
    }
    if (c == "22"){
      cluster.expression.22[,n] = data.expression.rna[,"22"]
    }
    if (c == "23"){
      cluster.expression.23[,n] = data.expression.rna[,"23"]
    }
    if (c == "24"){
      cluster.expression.24[,n] = data.expression.rna[,"24"]
    }
    if (c == "25"){
      cluster.expression.25[,n] = data.expression.rna[,"25"]
    }
    if (c == "26"){
      cluster.expression.26[,n] = data.expression.rna[,"26"]
    }
    if (c == "27"){
      cluster.expression.27[,n] = data.expression.rna[,"27"]
    }
    if (c == "28"){
      cluster.expression.28[,n] = data.expression.rna[,"28"]
    }
    if (c == "29"){
      cluster.expression.29[,n] = data.expression.rna[,"29"]
    }
    if (c == "30"){
      cluster.expression.30[,n] = data.expression.rna[,"30"]
    }
    }
}
              
write.table(cluster.expression.0, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster0.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.1, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster1.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.2, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster2.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.3, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster3.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.4, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster4.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.5, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster5.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.6, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster6.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.7, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster7.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.8, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster8.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.9, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster9.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.10, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster10.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.11, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster11.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.12, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster12.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.13, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster13.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.14, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster14.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.15, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster15.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.16, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster16.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.17, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster17.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.18, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster18.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.19, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster19.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.20, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster20.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.21, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster21.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.22, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster22.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.23, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster23.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.24, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster24.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.25, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster25.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.26, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster26.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.27, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster27.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.28, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster28.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.29, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster29.txt", row.names=TRUE, col.names=TRUE, sep='\t')
write.table(cluster.expression.30, "/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster30.txt", row.names=TRUE, col.names=TRUE, sep='\t')


