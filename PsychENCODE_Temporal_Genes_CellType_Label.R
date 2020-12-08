
brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 2.tsv"), header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/temporal_marker_genes.txt"), header=TRUE, row.names=1, sep="\t")

head(brain_genes)

print(brain_genes[1,"Cell.Type"])

celltypes <- character()
index = 1
for (gene in diff_expressed$name) {
  if (gene %in% brain_genes$Human.Gene) {
    print(gene)
    celltypes <- c(celltypes, brain_genes[index,"Cell.Type"])
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
  index = index + 1
}

diff_expressed$cell_type <- celltypes

write.table(diff_expressed, file = path.expand("~/temporal_marker_genes_celltypes.txt"), sep = '\t')
