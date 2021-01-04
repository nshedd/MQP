library(ArchR)
library(Seurat)

time = "PFC_Dev_Analysis"

proj <- loadArchRProject(path = time)

path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

seRNA.data = read.table(path1, header=TRUE, row.names=1)

seRNA <- CreateSeuratObject(counts = seRNA.data, project = "seRNA3k", min.cells = 3, min.features = 200)

#names(cluster_letters)=rownames(seRNA@meta.data)
cluster_letters <- Idents(object = seRNA)
names(cluster_letters) <- colnames(x = seRNA)
seRNA <- AddMetaData(
  object = seRNA,
  metadata = cluster_letters,
  col.name = 'letter.idents'
)
print("Added metadata")

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  groupRNA = "letter.idents",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)

cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))

unique(unique(proj$predictedGroup_Un))

write.table(cM, file = path.expand("~/git/MQP/PFC_Dev_Analysis/Plots/ConfusionMatrix-PFC-Integrated.txt"), sep = '\t')

#pal <- paletteDiscrete(values = colData(seRNA)$letter.idents)

p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un",
)

plotPDF(p1, name = "UMAP-PFC-Integrated", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

labelOld <- rownames(cM)

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
print("new labels")
labelNew

proj$Clusters2 <- mapLabels(proj$Clusters, newLabels = labelNew, oldLabels = labelOld)

p2 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "Clusters2",
)

plotPDF(p2, name = "UMAP-PFC-CLusters-Renamed", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)
