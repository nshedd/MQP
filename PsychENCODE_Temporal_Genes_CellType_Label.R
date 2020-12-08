
brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 2.tsv"), header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/temporal_marker_genes.txt"), header=TRUE, row.names=1, sep="\t")

head(brain_genes)

celltypelist = brain_genes$Cell.type
print(celltypelist[10])

celltypes <- character()
index = 1
for (gene in diff_expressed$name) {
  if (gene %in% brain_genes$Human.Gene) {
    print(celltypelist[index])
    celltypes <- c(celltypes, celltypelist[index])
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
  index = index + 1
}

diff_expressed$cell_type <- celltypes

write.table(diff_expressed, file = path.expand("~/temporal_marker_genes_celltypes.txt"), sep = '\t')
