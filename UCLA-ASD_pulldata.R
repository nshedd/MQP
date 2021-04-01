library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

group = c('CTL','ASD')
sample = unique(BA4.6$Sample)
cluster = unique(Idents(BA4.6))

for (c in cluster){
  for (s in sample) {
    data <- as.matrix(GetAssayData(BA4.6, slot = "counts"))
    data <- subset(data, idents = c)
    data <- subset(data, subset = "Sample"==s)
    head(data)
  }
}

