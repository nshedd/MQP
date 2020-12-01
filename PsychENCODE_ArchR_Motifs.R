library(ArchR)

time = "Temporal_Dev_Analysis"

proj <- loadArchRProject(path = time)

#if("Motif" %ni% names(proj@peakAnnotation)){
#  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
#}

#saveArchRProject(ArchRProj = proj)

proj <-addBgdPeaks(proj, method="chromVAR")

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)

plotVarDev

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = proj, addDOC = FALSE)

saveArchRProject(ArchRProj = proj)
