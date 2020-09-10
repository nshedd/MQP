df = readRDS("/data/projects/psychencode/df.rds")
library(umap)
library(ggplot2)

plotumap <- function(data) {
  data_umap = umap(data)
  data_layout = data.frame(data_umap$layout)
  ggplot(data=data_layout, aes(x=X1, y=X2,)) +geom_point()
}

plotumap(df)
