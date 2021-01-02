library(ArchR)
library(Seurat)

time = "PFC_Dev_Analysis"

proj <- loadArchRProject(path = time)

path1 = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

seRNA.data = read.table(path1, header=TRUE, row.names=1)

seRNA <- CreateSeuratObject(counts = seRNA.data, project = "seRNA3k", min.cells = 3, min.features = 200)

#names(cluster_letters)=rownames(seRNA@meta.data)
cluster_letters <- LETTERS[Idents(object = seRNA)]
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

proj

cM <- as.matrix(confusionMatrix(proj$Clusters, proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))

unique(unique(proj$predictedGroup_Un))
