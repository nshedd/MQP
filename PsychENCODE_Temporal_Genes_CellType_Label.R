
brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 2.tsv"), header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/temporal_marker_genes.txt"), header=TRUE, row.names=1, sep="\t")

head(brain_genes)

celltypes <- character()
for (gene in diff_expressed$name) {
  if (gene %in% brain_genes$Human.Gene) {
    celltypes <- c(celltypes, brain_genes$Cell.Type)
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

write.table(diff_expressed, file = path.expand("~/temporal_marker_genes_celltypes.txt"), sep = '\t')
