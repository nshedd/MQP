df = readRDS("/data/projects/psychencode/df.rds")
library(umap)
library(ggplot2)

plotumap <- function(data) {
  data_umap = umap(data)
  data_layout = data.frame(data_umap$layout)
    
  png("UMAPplot.png")
  ggplot(data=data_layout, aes(x=X1, y=X2, color=types[,2])) +geom_point()
  dev.off
}

plotumap(df)
