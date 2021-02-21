library(ArchR)

addArchRGenome("hg38")

args <- commandArgs(trailingOnly = TRUE)
time = "Lake_Integration"
print(time)

bam = "/home/sheddn/LakeTHSseq/data/TGCAGCTA_TCTCTCCG_TTTACC.unique.R7P3P4H_noAlt.FC9.rm_mm.bam.txt"
key = "FCX1"

ArrowFiles = createArrowFiles( inputFiles = bam, sampleNames = key,
filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE, 
gsubExpression=":.*")

proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)
rdhss = import("/data/projects/encode/Registry/V2/GRCh38/GRCh38-rDHSs.bed")
