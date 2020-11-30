library(ArchR)

time = "Temporal_Dev_Analysis"

proj <- loadArchRProject(path = time)

if("Motif" %ni% names(projHeme5@peakAnnotation)){
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
}

proj <-addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

