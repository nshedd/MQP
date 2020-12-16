matrixa = readRDS(path.expand("~/GSE97930_VisualCortex.RDS"))
matrixb = readRDS(path.expand("~/GSE97930_FrontalCortex.RDS"))
#matrixc = readRDS(path.expand("~/GSE97930_CerebellarHem.RDS"))

combinedmatrix = merge(matrixa, matrixb, by=0)

saveRDS(combinedmatrix, file = path.expand("~/GSE97930_All.RDS"))
