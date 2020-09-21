
df = read.table("/data/zusers/pratth/sc/atac/PBMC/PMBC.matrix.txt", sep="\t", header=TRUE)

library(umap)
library(ggplot2)

tdf = t(df)

plotumap <- function(data) {
  data_umap = umap(data, min_dist=0.02, n_neighbors=5)
  data_layout = data.frame(data_umap$layout)
  
  ggplot(data=data_layout, aes(x=X1, y=X2)) +geom_point()
  ggsave("ImmuneUMAP.png", device=)
}

plotumap(tdf)
