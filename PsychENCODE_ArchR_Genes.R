library(ArchR)

time = "Temporal_Dev_Analysis"

proj <- loadArchRProject(path = time)

proj <- addCoAccessibility(
  ArchRProj = proj,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = proj,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)

metadata(cA)[[1]]
