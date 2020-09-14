
df = readRDS("/data/projects/psychencode/df.rds")
types = read.table("cell_types.tsv", sep ="\t", header=FALSE)
 
library(umap)
library(ggplot2)

dl = types[,2]

tdf = t(df)

plotumap <- function(data, labels) {
  data_umap = umap(data)
  data_layout = data.frame(data_umap$layout)
  
  ggplot(data=data_layout, aes(x=X1, y=X2, color=dl)) +geom_point()
  ggsave("UMAPplot.png", device=)
}

plotumap(tdf, types)
