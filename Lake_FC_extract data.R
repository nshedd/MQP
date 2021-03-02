library(dplyr)
library(Seurat)

path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "Datavis", min.cells = 3, min.features = 200)

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)

FrontalCortex <- FindVariableFeatures(object = FrontalCortex)

all_cells <- rownames(FrontalCortex)
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))

#FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
#FrontalCortex <- FindClusters(FrontalCortex, resolution = 1)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20)

saveRDS(FrontalCortex,"~/Lake/FrontalCortex/Datavis_project.RDS")

expression = GetAssayData(object = FrontalCortex, slot = "scale.data")
write.table(expression, file = path.expand("~/Lake/FrontalCortex/datavis-expressionmatrix.txt"), sep="\t")

embeddings = as.data.frame(pbmc[["umap"]]@cell.embeddings)
write.table(embeddings, file = path.expand("~/Lake/FrontalCortex/datavis-embeddings.txt"), sep="\t")
