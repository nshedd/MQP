library(ArchR)
library(dplyr)
library(Seurat)
library(patchwork)

sample_names = scan(what = path.expand("~/sample-names.txt"), sep = "\t")
sample_names

addArchRGenome("hg19")
names("/data/zusers/pratth/sc/PEC/Temporal_GW20.bam")

ArrowFiles <- createArrowFiles(
  inputFiles = "/data/zusers/pratth/sc/PEC/Temporal_GW20.bam",
  sampleNames = names("/data/zusers/pratth/sc/PEC/Temporal_GW20.bam"),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

print("Done loading Archr file")

ArrowFiles

psychENCODE_whole <- CreateSeuratObject(counts = ArrowFiles, project = "psychENCODE_whole", min.cells = 3, min.features = 200)
psychENCODE_whole
