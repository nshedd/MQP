brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 2.tsv"), header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/GSE97930_CerebellarHem_differentiallyexpressed.txt"), header=TRUE, row.names=1, sep="\t")


celltypes <- character()
for (gene in diff_expressed$gene) {
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

write.table(diff_expressed, file = path.expand("~/GSE97930_CerebellarHem_differentiallyexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/GSE97930_CerebellarHem_differentiallyexpressed_condensed.txt"), sep = '\t')