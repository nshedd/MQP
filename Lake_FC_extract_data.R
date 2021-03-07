library(dplyr)
library(Seurat)
library(ggplot2)

path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

FrontalCortex <- CreateSeuratObject(counts = matrix, project = "Datavis", min.cells = 3, min.features = 200)

FrontalCortex[["percent.mt"]] <- PercentageFeatureSet(FrontalCortex, pattern = "^MT-")

FrontalCortex <- NormalizeData(FrontalCortex, normalization.method = "LogNormalize", scale.factor = 10000)

FrontalCortex <- FindVariableFeatures(object = FrontalCortex)

all_cells <- rownames(FrontalCortex)
FrontalCortex <- ScaleData(FrontalCortex, features = all_cells)

FrontalCortex <- RunPCA(FrontalCortex, features = VariableFeatures(object = FrontalCortex))

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 1)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20, metric="euclidean")

plot = DimPlot(FrontalCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/FrontalCortex_umap.png"), device=)

embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)

write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/FrontalCortex_embeddings.txt"), sep=",")

expression = GetAssayData(object = FrontalCortex, slot = "scale.data")
print("transposing...")
expression = t(expression)
print("done with transpose")
head(expression)
write.table(expression, file = path.expand("/data/rusers/sheddn/datavis4/FrontalCortex_expression.txt"), sep=",")
print("done writing table first time")

# markergenes = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B", "GAD1", "GAD2", "SLC6A1", "CLDN11", "MOG", "MOBP", "MBP", "SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4",
#                 "COBLL1", "DUSP1", "FLT1", "PCDH15", "OLIG1", "PCDH15", "OLIG1", "APBB1IP", "P2RY12", "APBB1IP", "P2RY12", "RYR1", "RELN", "GRM4", "RBFOX3")
# 
# expression = expression[markergenes]
write.table(expression, file = path.expand("/data/rusers/sheddn/datavis4/FrontalCortex_expression.txt"), sep=",")




path1 = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

matrix = read.table(path1, header=TRUE, row.names=1)

VisualCortex <- CreateSeuratObject(counts = matrix, project = "Datavis", min.cells = 3, min.features = 200)

VisualCortex[["percent.mt"]] <- PercentageFeatureSet(VisualCortex, pattern = "^MT-")

VisualCortex <- NormalizeData(VisualCortex, normalization.method = "LogNormalize", scale.factor = 10000)

VisualCortex <- FindVariableFeatures(object = VisualCortex)

all_cells <- rownames(VisualCortex)
VisualCortex <- ScaleData(VisualCortex, features = all_cells)

VisualCortex <- RunPCA(VisualCortex, features = VariableFeatures(object = VisualCortex))

VisualCortex <- FindNeighbors(VisualCortex, dims = 1:20)
VisualCortex <- FindClusters(VisualCortex, resolution = 1)

VisualCortex <- RunUMAP(VisualCortex, dims = 1:20, metric="euclidean")

plot = DimPlot(VisualCortex, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/VisualCortex_umap.png"), device=)

embeddings = as.data.frame(VisualCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(VisualCortex)

write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/VisualCortex_embeddings.txt"), sep=",")

expression = GetAssayData(object = VisualCortex, slot = "scale.data")
print("transposing...")
expression = t(expression)
print("done with transpose")

# markergenes = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B", "GAD1", "GAD2", "SLC6A1", "CLDN11", "MOG", "MOBP", "MBP", "SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4",
#                 "COBLL1", "DUSP1", "FLT1", "PCDH15", "OLIG1", "PCDH15", "OLIG1", "APBB1IP", "P2RY12", "APBB1IP", "P2RY12", "RYR1", "RELN", "GRM4", "RBFOX3")
# 
# expression = expression[markergenes]
write.table(expression, file = path.expand("/data/rusers/sheddn/datavis4/VisualCortex_expression.txt"), sep=",")



path1 = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrix = read.table(path1, header=TRUE, row.names=1)

CerebellarHem <- CreateSeuratObject(counts = matrix, project = "Datavis", min.cells = 3, min.features = 200)

CerebellarHem[["percent.mt"]] <- PercentageFeatureSet(CerebellarHem, pattern = "^MT-")

CerebellarHem <- NormalizeData(CerebellarHem, normalization.method = "LogNormalize", scale.factor = 10000)

CerebellarHem <- FindVariableFeatures(object = CerebellarHem)

all_cells <- rownames(CerebellarHem)
CerebellarHem <- ScaleData(CerebellarHem, features = all_cells)

CerebellarHem <- RunPCA(CerebellarHem, features = VariableFeatures(object = CerebellarHem))

CerebellarHem <- FindNeighbors(CerebellarHem, dims = 1:20)
CerebellarHem <- FindClusters(CerebellarHem, resolution = 1)

CerebellarHem <- RunUMAP(CerebellarHem, dims = 1:20, metric="euclidean")

plot = DimPlot(CerebellarHem, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
ggsave(path.expand("/data/rusers/sheddn/datavis4/CerebellarHem_umap.png"), device=)

embeddings = as.data.frame(CerebellarHem[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(CerebellarHem)

write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/CerebellarHem_embeddings.txt"), sep=",")

expression = GetAssayData(object = CerebellarHem, slot = "scale.data")
print("transposing...")
expression = t(expression)
print("done with transpose")

# markergenes = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B", "GAD1", "GAD2", "SLC6A1", "CLDN11", "MOG", "MOBP", "MBP", "SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4",
#                 "COBLL1", "DUSP1", "FLT1", "PCDH15", "OLIG1", "PCDH15", "OLIG1", "APBB1IP", "P2RY12", "APBB1IP", "P2RY12", "RYR1", "RELN", "GRM4", "RBFOX3")
# 
# expression = expression[markergenes]
write.table(expression, file = path.expand("/data/rusers/sheddn/datavis4/CerebellarHem_expression.txt"), sep=",")
