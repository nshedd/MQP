library(ArchR)

#GRCh38
addArchRGenome("hg19")

args <- commandArgs(trailingOnly = TRUE)
time = args[1]
key = c("Temporal_GW20")
bam = c("/data/zusers/pratth/sc/PEC/Temporal_GW20.tsv.gz")

ArrowFiles = createArrowFiles( inputFiles = bam, sampleNames = key,
                               filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE, 
                               #gsubExpression=":.*"
                               )
ArrowFiles


proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)
