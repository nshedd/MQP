library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

group = c('CTL','ASD')
sample = c('BW11','BW12','BW13','BW15','BW2','BW23','BW25','BW30','BW35','BW46','BW5','BW6','BW8','BW1','BW10','BW14','BW16','BW17','BW22','BW24','BW26','BW27','BW3','BW31','BW32','BW4','BW59','BW60','BW7','BW9')
cluster = range(0,38)

for (c in cluster){
  for (s in sample) {
    data <- as.matrix(GetAssayData(BA4.6, slot = "counts")[, WhichCells(BA4.6, ident=c)])
    data <- subset(data, subset = "Sample"==s)
  }
}

head(data)
