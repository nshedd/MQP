library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

install.packages('BiocManager')
BiocManager::install('limma') 

set2umap = readRDS(file=path.expand("~/GSM3722075_PBMC_Rep3_fragments_Seurat.rds"))

T_cells = read.table("/data/zusers/pratth/ATAC/specific-elements/top-10k/unstimulated_T-cells.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
B_cells = read.table("/data/zusers/pratth/ATAC/specific-elements/top-10k/B-cell.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
M_cells = read.table("/data/zusers/pratth/ATAC/specific-elements/top-10k/myeloid-cells.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
print("read tables")

set2umap.markers <- FindAllMarkers(set2umap, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
diff_expressed = set2umap.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

write.table(diff_expressed, file = path.expand("~/set2_differentiallyexpressed.txt"), sep="\t")
print("wrote table once")

celltypes <- character()
for (marker in diff_expressed$gene) {
  if (marker %in% T_cells) {
    celltypes <- c(celltypes, "t-cell")
  } else if (marker %in% B_cells) {
    celltypes <- c(celltypes, "b-cell")
  } else if (marker %in% M_cells) {
    celltypes <- c(celltypes, "m-cell")
  } else {
    celltypes <- c(celltypes, "unknown")
  }
}

diff_expressed$cell_type <- celltypes

print("finished matching cell types")

write.table(diff_expressed, file = path.expand("~/set2_differentiallyexpressed.txt"), sep="\t")
print("wrote table twice")
