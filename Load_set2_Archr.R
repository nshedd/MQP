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

p=getEmbedding(ArchRProj = proj, embedding = "UMAP_Peak")
t=getEmbedding(ArchRProj = proj, embedding = "UMAP_Tile")

png(paste0("UMAP-Peak-",time,".png"), width=10, height=8, units="in", res=300)
ggplot(p, aes(x=p[,1], y=p[,2], color=proj@cellColData$Clusters_Peak))+
  geom_point()#+scale_color_manual(values=c("red4","limegreen","green4","red3","darkorchid","violetred2","royalblue"))
dev.off()

png(paste0("UMAP-Tile-",time,".png"), width=10, height=8, units="in", res=300)
ggplot(t, aes(x=t[,1], y=t[,2], color=proj@cellColData$Clusters_Tile))+
  geom_point()#+scale_color_manual(values=c("red4","limegreen","green4","darkorchid","violetred2","royalblue"))
dev.off()

saveArchRProject(ArchRProj = proj, outputDirectory = time, overwrite = TRUE, load = TRUE)
