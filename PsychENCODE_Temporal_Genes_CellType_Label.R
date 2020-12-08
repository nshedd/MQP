
brain_genes = read.table("Zlab single-cell marker genes - Brain.tsv 2", header=TRUE, sep="\t")

diff_expressed = read.table(path.expand("~/temporal_marker_genes.txt"), header=TRUE, row.names=1, sep="\t")

print(

celltypes <- character()
for (gene in diff_expressed$name) {
  if (gene %in% brain_genes$Human_Gene) {
    celltypes <- c(celltypes, brain_genes$Human_Gene)
  }
  else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

write.table(diff_expressed, file = path.expand("~/temporal_marker_genes_celltypes.txt"), sep = '\t')
