library(ArchR)

library(ArchR)

time = "Temporal_Dev_Analysis"

proj <- loadArchRProject(path = time)

if("Motif" %ni% names(proj@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = proj,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)