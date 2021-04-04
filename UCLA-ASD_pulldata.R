library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

data.expression = AverageExpression(BA4.6, slot="counts", group.by=Sample)
head(data.expression[,1:5])
q()


group = c('CTL','ASD')
sample = unique(BA4.6$Sample)
cluster = unique(Idents(BA4.6))
               
for (c in cluster) {
  data.expression = AverageExpression(data, slot="counts", group.by=Sample)
  # data.expresssion = as.matrix(data.expression[["RNA"]])
}
              
head(data.expression[,1:5])
