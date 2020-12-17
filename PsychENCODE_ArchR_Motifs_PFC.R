library(ArchR)

time = "PFC_Dev_Analysis"

proj <- loadArchRProject(path = time)

if("Motif_ENCODE" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "encode", name = "Motif_ENCODE")
}

markersPeaks <- getMarkerFeatures(ArchRProj = proj, useMatrix = "PeakMatrix", groupBy = "Clusters",
                                  bias = c("TSSEnrichment", "log10(nFrags)"),testMethod = "wilcoxon")

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif_ENCODE",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap_ENCODE_PFC", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)
