iris = read.table("iris.csv", header=TRUE, sep=",")

irisData = iris[,1:4]
irisLabels = iris[,5]

library(umap)
library(ggplot2)

plotumap <- function(data, labels) {
  data_umap = umap(data)
  data_layout = data.frame(data_umap$layout)
  
  png("UMAPplot.png")
  ggplot(data=data_layout, aes(x=X1, y=X2, color=labels)) +geom_point()
  dev.off()
}

plotumap(irisData, irisLabels)
