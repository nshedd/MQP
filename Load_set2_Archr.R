library(ArchR)

addArchRGenome("hg38")

args <- commandArgs(trailingOnly = TRUE)
time = args[1]
bam = "/data/zusers/pratth/sc/PEC/Temporal_GW20.sorted.fragments.tsv.gz"
key = "temporal"

ArrowFiles = createArrowFiles( inputFiles = bam, sampleNames = key,
  filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE, 
  gsubExpression=":.*")

proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)
rdhss = import("/data/projects/encode/Registry/V2/GRCh38/GRCh38-rDHSs.bed")

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

proj <- addClusters(input = proj, reducedDims = "IterativeLSI")

proj <- addUMAP(ArchRProj = proj, reducedDims = "IterativeLSI")

p <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")

plotPDF(p, name = "Plot-UMAP-Sample-Clusters.pdf",
        ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = proj, outputDirectory = time, overwrite = TRUE, load = TRUE)
