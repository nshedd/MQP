brain_genes = read.table(path.expand("~/Zlab single-cell marker genes - Brain 2.tsv"), header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/GSE97930_visualcortex_differentiallyexpressed.txt"), header=TRUE, row.names=1, sep="\t")


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

diff_expressed_condensed_pres = diff_expressed_condensed[c("cluster", "gene", "cell_type")]

write.table(diff_expressed, file = path.expand("~/GSE97930_visualcortex_differentiallyexpressed.txt"), sep = '\t')

write.table(diff_expressed_condensed, file = path.expand("~/GSE97930_visualcortex_differentiallyexpressed_condensed.txt"), sep = '\t')

write.table(diff_expressed_condensed_pres, file = path.expand("~/GSE97930_visualcortex_differentiallyexpressed_condensed_pres.txt"), sep = '\t')
