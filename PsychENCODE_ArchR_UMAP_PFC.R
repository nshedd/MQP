library(ArchR)

addArchRGenome("hg38")

time = "PFC_Dev_Analysis"
print(time)
bam = "/data/zusers/pratth/sc/PEC/Twin31_PFC.bam.tsv.tsv.gz"
key = "prefrontal"

ArrowFiles = createArrowFiles( inputFiles = bam, sampleNames = key,
  filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE, 
  gsubExpression=":.*")

proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)
rdhss = import("/data/projects/encode/Registry/V2/GRCh38/GRCh38-rDHSs.bed")

proj = addPeakSet(ArchRProj = proj, peakSet = rdhss, force = FALSE)
proj = addPeakMatrix(proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "PeakMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, nNeighbors = 10, reducedDims = "IterativeLSI")

p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")		

plotPDF(p, name = "Plot-UMAP-PFC-Clusters.pdf",		
         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)		

markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
write.table(markerList, file = path.expand("~/git/MQP/PFC_Dev_Analysis/Tables/pfc_marker_peaks.txt"), sep = '\t')

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-PFC", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)
