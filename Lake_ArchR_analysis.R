library(ArchR)

addArchRGenome("hg38")

args <- commandArgs(trailingOnly = TRUE)
time = "/data/rusers/sheddn/Lake_ths-seq/Lake_Integration_Analysis"
print(time)
fragment = "/data/rusers/sheddn/Lake_ths-seq/fragments.lake.no-underscores.tsv.gz"
key = "Lake"

reformatFragmentFiles(
  fragmentFiles = fragment,
  checkChrPrefix = getArchRChrPrefix()
)

ArrowFiles = createArrowFiles( inputFiles = fragment, sampleNames = key,
  filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE, force=TRUE,
  gsubExpression=":.*")

proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)
print(proj$cellNames)
rdhss = import("/data/projects/encode/Registry/V2/GRCh38/GRCh38-rDHSs.bed")
proj = addPeakSet(ArchRProj = proj, peakSet = rdhss, force = TRUE)
proj = addPeakMatrix(proj)

print(proj$cellNames)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)	

proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")		

plotPDF(p, name = "Plot-Lake-UMAP-Clusters.pdf",		
         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)		

markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
write.table(markerList, file = "/data/rusers/sheddn/Lake_ths-seq/overenriched_peaks.txt", sep = '\t')

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Heatmap-Lake", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)


markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

write.table(markerList, file ="/data/rusers/sheddn/Lake_ths-seq/overexpressed_genes.txt", sep = '\t')

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

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

write.table(cM, file = path.expand("/data/rusers/sheddn/Lake_ths-seq/ConfusionMatrix-PFC-Integrated.txt"), sep = '\t')

#pal <- paletteDiscrete(values = colData(seRNA)$letter.idents)

p1 <- plotEmbedding(
    proj, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un",
)

plotPDF(p1, name = "UMAP-PFC-Integrated", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)
