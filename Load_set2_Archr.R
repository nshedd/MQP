library(ArchR)

sample_names = scan(what = path.expand("~/sample-names.txt"), sep = "\t")
print(sample_names)

ArrowFiles <- createArrowFiles(
  inputFiles = "/data/zusers/pratth/sc/atac/GSM3722075_PBMC_Rep3_fragments.tsv.gz.rDHS.matrix.tsv",
  sampleNames = sample_names,
  #filterTSS = 4, #Dont set this too high because you can always increase later
  #filterFrags = 1000, 
  #addTileMat = TRUE,
  #addGeneScoreMat = TRUE
)

print("Done loading Archr file")
