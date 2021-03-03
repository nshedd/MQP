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

FrontalCortex <- FindNeighbors(FrontalCortex, dims = 1:20)
FrontalCortex <- FindClusters(FrontalCortex, resolution = 1)

FrontalCortex <- RunUMAP(FrontalCortex, dims = 1:20)

embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
head(embeddings)
embeddings$celltype = Idents(FrontalCortex)

embeddings = as.data.frame(FrontalCortex[["umap"]]@cell.embeddings)
write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings.txt"), sep="\t")

avg_expression_full = AverageExpression(FrontalCortex)[["RNA"]]

head(avg_expression_full)

avg_expression_full = as.data.frame(avg_expression_full)

num_clusters = nrow(avg_expression_full)

ex = (avg_expression_full$SLC17A7 + avg_expression_full$GRIN1 + avg_expression_full$GRIN2B + avg_expression_full$SATB2) / 4
print(avg_expression_full$SLC17A7)
inh = (avg_expression_full$GAD1 + avg_expression_full$GAD2 + avg_expression_full$SLC6A1) / 3
print(inh)
oli = (avg_expression_full$CLDN11 + avg_expression_full$MOG + avg_expression_full$MOBP + avg_expression_full$MBP) / 4
print(oli)
ast = (avg_expression_full$SLC1A2 + avg_expression_full$SLC1A3 + avg_expression_full$LC4A4 + avg_expression_full$GLUL + avg_expression_full$AQP4) / 5
print(ast)
mic = (avg_expression_full$APBB1IP + avg_expression_full$P2RY12) / 2
print(mic)
opc = (avg_expression_full$PCDH15 + avg_expression_full$OLIG1) / 2
print(opc)
end = (avg_expression_full$COBLL1 + avg_expression_full$DUSP1 + avg_expression_full$FLT1) / 3
print(end)

avg_expression = data.frame('cluster'= 1:num_clusters,
                            'Excitatory neuron'= ex,
                            'Inhibitory neuron'= inh,
                            'Oligodendrocyte'= oli,
                            'Astrocyte'= ast,
                            'Microglia'= mic,
                            'Oligodendrocyte precursor'= opc,
                            'Endothelial cell'= end)
                            
write.table(expression, file = path.expand("/data/rusers/sheddn/datavis4/expressionmatrix.txt"), sep="\t")

# expression = GetAssayData(object = FrontalCortex, slot = "scale.data")
# 
# marker_genes = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B", "GAD1", "GAD2", "SLC6A1", "CLDN11", "MOG", "MOBP", "MBP", "SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4"
#                 "COBLL1", "DUSP1", "FLT1", "PCDH15", "OLIG1", "PCDH15", "OLIG1", "APBB1IP", "P2RY12", "APBB1IP", "P2RY12", "RYR1", "RELN", "GRM4", "RBFOX3")
# 
# expression = expression[markergenes,]
# write.table(expression, file = path.expand("~/Lake/FrontalCortex/datavis-expressionmatrix.txt"), sep="\t")

