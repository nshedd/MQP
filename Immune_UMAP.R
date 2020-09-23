
df = read.table("/data/zusers/pratth/sc/atac/PBMC/PBMC.matrix.txt", sep="\t", header=TRUE, row.names=1)

library(umap)
library(ggplot2)

tdf = t(df)

plotumap <- function(data) {
  data_umap = umap(data)
  data_layout = data.frame(data_umap$layout)
  
  ggplot(data=data_layout, aes(x=X1, y=X2)) +geom_point()
  ggsave("ImmuneUMAP.png", device=)
}

plotumap(df)
