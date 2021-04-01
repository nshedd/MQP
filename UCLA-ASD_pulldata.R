library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

group = c('CTL','ASD')
sample = unique(BA4.6$Sample)
cluster = unique(Idents(BA4.6))

for (s in sample) {
  data <- subset(BA4.6, subset = Sample==s)
  data.expression = AverageExpression(data, slot="counts")
  data.expresssion = as.matrix(data.expression[["RNA"]])
}

head(data.expression)
