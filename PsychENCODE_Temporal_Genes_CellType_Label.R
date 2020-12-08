
brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 2.tsv"), header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/temporal_marker_genes.txt"), header=TRUE, row.names=1, sep="\t")

celltypelist = brain_genes$Cell.type
print(celltypelist[277])

celltypes <- character()
ind = 1
for (gene in diff_expressed$name) {
  if (gene %in% brain_genes$Human.Gene) {
    print(ind)
    print(celltypelist[ind])
    celltypes <- c(celltypes, celltypelist[ind])
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
  ind = ind + 1
}

diff_expressed$cell_type <- celltypes

write.table(diff_expressed, file = path.expand("~/temporal_marker_genes_celltypes.txt"), sep = '\t')
