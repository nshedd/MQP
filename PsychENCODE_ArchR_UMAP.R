library(ArchR)

addArchRGenome("hg38")

args <- commandArgs(trailingOnly = TRUE)
time = args[1]
print(time)
bam = "/data/zusers/pratth/sc/PEC/Temporal_GW20.sorted.fragments.tsv.gz"
key = "temporal"

ArrowFiles = createArrowFiles( inputFiles = bam, sampleNames = key,
  filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE, 
  gsubExpression=":.*")

proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)
rdhss = import("/data/projects/encode/Registry/V2/GRCh38/GRCh38-rDHSs.bed")

proj = addPeakSet(ArchRProj = proj, peakSet = rdhss, force = FALSE)
proj = addPeakMatrix(proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")		

ggAlignPlots(p1, p2, type = "h")		

plotPDF(p1, p2, name = "Plot-UMAP-Temporal-Clusters-2.pdf",		
         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)		

markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
write.table(markerList, file = path.expand("~/temporal_marker_peaks.txt"), sep = '\t')

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap-Temporal", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)
