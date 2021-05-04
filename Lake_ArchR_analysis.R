library(ArchR)

addArchRGenome("hg38")

#marker_gene_table = read.table("~/Zlab single-cell marker genes - Brain 4.tsv", header=TRUE, sep="\t")
#head(marker_gene_table)

args <- commandArgs(trailingOnly = TRUE)
time = "/data/rusers/sheddn/Lake_ths-seq/Lake_Integration_Analysis"
print(time)
fragment = "/data/rusers/sheddn/Lake_ths-seq/Lake-040521.tsv.gz"
key = "Lake"

reformatFragmentFiles(
  fragmentFiles = fragment,
  checkChrPrefix = getArchRChrPrefix()
)

ArrowFiles = createArrowFiles( inputFiles = fragment, sampleNames = key,
  minTSS = 2, minFrags = 100, addTileMat = TRUE, addGeneScoreMat = TRUE,
  gsubExpression=":.*")

proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))

p <- ggPoint(
    x = df[,1], 
    y = df[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projHeme1, addDOC = FALSE)

q()

print(proj$cellNames)
rdhss = import("/data/projects/encode/Registry/V2/GRCh38/GRCh38-rDHSs.bed")
proj = addPeakSet(ArchRProj = proj, peakSet = rdhss)
proj = addPeakMatrix(proj)

print(proj$cellNames)

proj <- addDoubletScores(
    input = proj,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)	

# proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "PeakMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI", method = "Seurat")

proj <- addUMAP(ArchRProj = proj, nNeighbors = 10, reducedDims = "IterativeLSI")

p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")		

plotPDF(p, name = "Plot-Lake-UMAP-Clusters.pdf",		
         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)		

# proj <- subsetArchRProject(proj)
# 
# saveArchRProject(ArchRProj = proj, overwrite = TRUE, load = TRUE)

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

markerList <- read.table("/data/rusers/sheddn/Lake_ths-seq/overexpressed_genes.txt", sep = '\t')

all_known_marker_genes = marker_gene_table$Human.Gene

intersection = intersect(markerList$name, all_known_marker_genes)

markerGenes  <- intersection

print(markerGenes)

write.table(markerGenes, file ="/data/rusers/sheddn/Lake_ths-seq/intesection.txt", sep = ',')

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

path1 = path.expand("/data/rusers/sheddn/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")

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

plotPDF(p1, name = "UMAP-Lake-Integrated", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

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

plotPDF(p2, name = "UMAP-Lake-Clusters-Integrated", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)
