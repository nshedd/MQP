
df = read.table("/data/zusers/pratth/sc/atac/PBMC", sep="\t", header=TRUE)

library(umap)
library(ggplot2)

tdf = t(df)

plotumap <- function(data, labels) {
  data_umap = umap(data, min_dist=0.02, n_neighbors=3)
  data_layout = data.frame(data_umap$layout)
  
  ggplot(data=data_layout, aes(x=X1, y=X2, color=labels)) +geom_point()
  ggsave("UMAPplot.png", device=)
}

plotumap(tdf, tdl)