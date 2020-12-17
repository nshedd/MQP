library(ArchR)

time = "PFC_Dev_Analysis"

proj <- loadArchRProject(path = time)

markersGS <- getMarkerFeatures(
    ArchRProj = proj, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

#write.table(markersGS, file = path.expand("~/pfc_marker_gene_score.txt"), sep = '\t')

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

write.table(markerList, file = path.expand("~/pfc_marker_genes.txt"), sep = '\t')

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap-PFC", width = 8, height = 6, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)
