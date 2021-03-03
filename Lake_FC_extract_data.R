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

write.table(embeddings, file = path.expand("/data/rusers/sheddn/datavis4/embeddings.txt"), sep=",")

avg_expression_full = t(AverageExpression(FrontalCortex)[["RNA"]])

write.table(avg_expression_full, file ="/data/rusers/sheddn/datavis4/expression_full.txt", sep="\t")
avg_expression_full = read.table("/data/rusers/sheddn/datavis4/expression_full.txt", header=TRUE, row.names=1)

num_clusters = nrow(avg_expression_full)

print(avg_expression_full$SLC17A7)

ex = rowMeans(avg_expression_full[,c("SLC17A7", "SATB2", "GRIN1", "GRIN2B")])
inh = rowMeans(avg_expression_full[,c("GAD1", "GAD2", "SLC6A1")])
oli = rowMeans(avg_expression_full[,c("CLDN11", "MOG", "MOBP", "MBP")])
ast = rowMeans(avg_expression_full[,c("SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4")])
mic = rowMeans(avg_expression_full[,c("APBB1IP", "P2RY12")])
opc = rowMeans(avg_expression_full[,c("PCDH15", "OLIG1")])
end = rowMeans(avg_expression_full[,c("COBLL1", "DUSP1", "FLT1")])
print(end)

avg_expression = data.frame('cluster'= c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,21),
                            'Excitatory neuron'= ex,
                            'Inhibitory neuron'= inh,
                            'Oligodendrocyte'= oli,
                            'Astrocyte'= ast,
                            'Microglia'= mic,
                            'Oligodendrocyte precursor'= opc,
                            'Endothelial cell'= end)

print(avg_expression)
                            
write.table(avg_expression, file = path.expand("/data/rusers/sheddn/datavis4/expressionmatrix.txt"), sep=",")

# expression = GetAssayData(object = FrontalCortex, slot = "scale.data")
# 
# marker_genes = c("SLC17A7", "SATB2", "GRIN1", "GRIN2B", "GAD1", "GAD2", "SLC6A1", "CLDN11", "MOG", "MOBP", "MBP", "SLC1A2", "SLC1A3", "SLC4A4", "GLUL", "AQP4",
#                 "COBLL1", "DUSP1", "FLT1", "PCDH15", "OLIG1", "PCDH15", "OLIG1", "APBB1IP", "P2RY12", "APBB1IP", "P2RY12", "RYR1", "RELN", "GRM4", "RBFOX3")
# 
# expression = expression[markergenes,]
# write.table(expression, file = path.expand("~/Lake/FrontalCortex/datavis-expressionmatrix.txt"), sep="\t")

