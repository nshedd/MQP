library(ArchR)

time = "Temporal_Dev_Analysis"

proj <- loadArchRProject(path = time)

proj <- addImputeWeights(proj)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "GeneScoreMatrix",
  embedding = "UMAP",
  imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))


plotPDF(plotList = p, 
        name = "Plot-UMAP-Temporal-Marker-Genes.pdf", 
        ArchRProj = proj, 
        addDOC = FALSE, width = 5, height = 5)
