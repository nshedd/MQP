library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

data.expression = AverageExpression(BA4.6, slot="counts")
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
cluster.expression.31 = data.expression.rna[,"31", drop=FALSE]
cluster.expression.32 = data.expression.rna[,"32", drop=FALSE]
cluster.expression.33 = data.expression.rna[,"33", drop=FALSE]
cluster.expression.34 = data.expression.rna[,"34", drop=FALSE]
cluster.expression.35 = data.expression.rna[,"35", drop=FALSE]
cluster.expression.36 = data.expression.rna[,"36", drop=FALSE]
cluster.expression.37 = data.expression.rna[,"37", drop=FALSE]
cluster.expression.38 = data.expression.rna[,"38", drop=FALSE]

write.table(data.expression.rna, "/data/rusers/sheddn/UCLA-ASD/data/BA4.6_averageexpression.txt")



group = c('CTL','ASD')
sample = unique(BA4.6$Sample)
cluster = unique(Idents(BA4.6))
               
for (s in sample) {
  data <- subset(BA4.6, subset= Sample==s)
  data.expression = AverageExpression(data, slot="counts")
  data.expression.rna = data.expression[["RNA"]]
  cluster.expression.0[,s] = data.expression.rna[,"0"]
  cluster.expression.1[,s] = data.expression.rna[,"1"]
  cluster.expression.2[,s] = data.expression.rna[,"2"]
  cluster.expression.3[,s] = data.expression.rna[,"3"]
  cluster.expression.4[,s] = data.expression.rna[,"4"]
  cluster.expression.5[,s] = data.expression.rna[,"5"]
  cluster.expression.6[,s] = data.expression.rna[,"6"]
  cluster.expression.7[,s] = data.expression.rna[,"7"]
  cluster.expression.8[,s] = data.expression.rna[,"8"]
  cluster.expression.9[,s] = data.expression.rna[,"9"]
  cluster.expression.10[,s] = data.expression.rna[,"10"]
  cluster.expression.11[,s] = data.expression.rna[,"11"]
  cluster.expression.12[,s] = data.expression.rna[,"12"]
  cluster.expression.13[,s] = data.expression.rna[,"13"]
  cluster.expression.14[,s] = data.expression.rna[,"14"]
  cluster.expression.15[,s] = data.expression.rna[,"15"]
  cluster.expression.16[,s] = data.expression.rna[,"16"]
  cluster.expression.17[,s] = data.expression.rna[,"17"]
  cluster.expression.18[,s] = data.expression.rna[,"18"]
  cluster.expression.19[,s] = data.expression.rna[,"19"]
  cluster.expression.20[,s] = data.expression.rna[,"20"]
  cluster.expression.21[,s] = data.expression.rna[,"21"]
  cluster.expression.22[,s] = data.expression.rna[,"22"]
  cluster.expression.23[,s] = data.expression.rna[,"23"]
  cluster.expression.24[,s] = data.expression.rna[,"24"]
  cluster.expression.25[,s] = data.expression.rna[,"25"]
  cluster.expression.26[,s] = data.expression.rna[,"26"]
  cluster.expression.27[,s] = data.expression.rna[,"27"]
  cluster.expression.28[,s] = data.expression.rna[,"28"]
  cluster.expression.29[,s] = data.expression.rna[,"29"]
  cluster.expression.30[,s] = data.expression.rna[,"30"]
  cluster.expression.31[,s] = data.expression.rna[,"31"]
  cluster.expression.32[,s] = data.expression.rna[,"32"]
  cluster.expression.33[,s] = data.expression.rna[,"33"]
  cluster.expression.34[,s] = data.expression.rna[,"34"]
  cluster.expression.35[,s] = data.expression.rna[,"35"]
  cluster.expression.36[,s] = data.expression.rna[,"36"]
  cluster.expression.37[,s] = data.expression.rna[,"37"]
  cluster.expression.38[,s] = data.expression.rna[,"38"]
}
              
head(cluster.expression.0)
