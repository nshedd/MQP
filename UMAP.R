

df = readRDS("/data/projects/psychencode/df.rds")
types = read.table("cell_types.tsv", sep ="\t", header=FALSE)
 
library(umap)
library(ggplot2)

ttypes = t(types)
dl = types[,2]
tdl = ttypes[2,]

tdf = t(df)

plotumap <- function(data, labels) {
  data_umap = umap(data, n_neighbors=3, min_dist = 0.05)
  data_layout = data.frame(data_umap$layout)
  
  ggplot(data=data_layout, aes(x=X1, y=X2, color=labels)) +geom_point()
  ggsave("UMAPplot_mindist05.png", device=)
}

plotumap(tdf, tdl)
