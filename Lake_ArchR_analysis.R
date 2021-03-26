library(ArchR)

addArchRGenome("hg38")

args <- commandArgs(trailingOnly = TRUE)
time = "/data/rusers/sheddn/Lake_ths-seq/Lake_Integration_Analysis"
print(time)
fragment = "/data/zusers/pratth/sc/fragments.lake.tsv.gz"
key = "Lake"

reformatFragmentFiles(
  fragmentFiles = fragment,
  checkChrPrefix = getArchRChrPrefix()
)

ArrowFiles = createArrowFiles( inputFiles = fragment, sampleNames = key,
  filterTSS = 4, filterFrags = 1000, addTileMat = TRUE, addGeneScoreMat = TRUE)

proj = ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = time, copyArrows = TRUE)
rdhss = import("/data/projects/encode/Registry/V2/GRCh38/GRCh38-rDHSs.bed")
