matrixa = readRDS(path.expand("~/GSE97930_VisualCortex.RDS"))
matrixb = readRDS(path.expand("~/GSE97930_FrontalCortex.RDS"))
matrixc = readRDS(path.expand("~/GSE97930_CerebellarHem.RDS"))
matrixab = readRDS(path.expand("~/GSE97930_All.RDS"))

combinedmatrix = merge(matrixa, matrixb, by=0)

row.names(combinedmatrix)<-combinedmatrix$Row.names
combinedmatrix[ , !("Row.names")]

combinedmatrix[0:10,0:10]

saveRDS(combinedmatrix, file = path.expand("~/GSE97930_All.RDS"))
