library(ArchR)

marker_gene_table = read.table("~/Zlab single-cell marker genes - Brain 4.tsv", header=TRUE, sep="\t")
head(marker_gene_table)

args <- commandArgs(trailingOnly = TRUE)
time = "/data/rusers/sheddn/Lake_ths-seq/Lake_Integration_Analysis"
print(time)
fragment = "/data/rusers/sheddn/Lake_ths-seq/Lake-040521.tsv.gz"
key = "Lake"

proj <- loadArchRProject(time)

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
