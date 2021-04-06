library(Seurat)
library(DESeq2)
## BA4/6

metaData = read.table("/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA4.6_metadata_cluster0.txt", sep='\t', header=TRUE)
head(metaData)

countData = read.table(file="/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA4.6_averageexpression_cluster0.txt", header=TRUE, sep="\t")
head(countData)

dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~1+Group)
dds <- DESeq(dds)

res <- results(dds)
head(results(dds, tidy=TRUE))

#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))

# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

ggsave("/data/rusers/sheddn/UCLA-ASD/plots/DEG_cluster0.png")

q()
for (i in range(1,38)) {
  metadata_file = paste("/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA4.6_metadata_cluster", i, ".txt", sep='')
  countData = read.table(file=metadata_file, sep="\t")
  
  expression_file = paste("/data/rusers/sheddn/UCLA-ASD/data/DEG_output/BA4.6_averageexpression_cluster", i, ".txt", sep='')
  countData = read.table(file=expression_file, sep="\t")
  head(countData)
  dds <- DESeqDataSetFromMatrix(countData=countData, 
                                colData=metaData, 
                                design=~1 + Sample + Group)
  dds <- DESeq(dds, full=~Group)
  
  res <- results(dds)
  head(results(dds, tidy=TRUE))
  
  #reset par
  par(mfrow=c(1,1))
  # Make a basic volcano plot
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-3,3)))
  
  # Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
  with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
  with(subset(res, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
  
  pic = paste("/data/rusers/sheddn/UCLA-ASD/plots/DEG_cluster", i, ".png", sep='')
  ggsave(pic)
}
