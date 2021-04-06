library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

group = c('CTL','ASD')
sample = unique(BA4.6$Sample)
cluster = unique(Idents(BA4.6))

for (c in cluster) {
  data <- subset(BA4.6, ident = c)
  expression.data <- GetAssayData(data, slot = "counts")
  expression_file <- paste("/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA4.6_averageexpression_cluster", c, ".txt", sep='')
  write.table(expression.data, expression_file, sep='\t')
  
  head(expression.data)
            
  df = data[[c("Group","Sample")]]
  meta_file <- paste("/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA4.6_metadata_cluster", c, ".txt", sep='')
  write.table(df, meta_file, sep='\t')
  
  head(df)
}


BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

sample = unique(BA9$Sample)
cluster = unique(Idents(BA9))

for (c in cluster) {
  data <- subset(BA9, ident = c)
  expression.data <- GetAssayData(data, slot = "counts")
  expression_file <- paste("/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_averageexpression_cluster", c, ".txt", sep='')
  write.table(expression.data, expression_file, sep='\t')
  
  head(expression.data)
            
  df = data[[c("Group","Sample")]]
  meta_file <- paste("/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA9_metadata_cluster", c, ".txt", sep='')
  write.table(df, meta_file, sep='\t')
  
  head(df)
}
