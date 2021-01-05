
brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 3.tsv"), header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/pfc_marker_genes.txt"), header=TRUE, row.names=1, sep="\t")

celltypelist = brain_genes$Cell.type
print(celltypelist[277])

celltypes <- character()
for (gene in diff_expressed$name) {
  if (gene %in% brain_genes$Human.Gene) {
    celltype = brain_genes$Cell.type[brain_genes$Human.Gene == gene]
    celltype = paste(celltype, collapse = ', ')
    celltypes <- c(celltypes, celltype)
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

diff_expressed_condensed = diff_expressed[diff_expressed$cell_type != "unknown",]

write.table(diff_expressed, file = path.expand("~/pfc_marker_genes_celltypes.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/pfc_marker_genes_celltypes_condensed.txt"), sep = '\t')
