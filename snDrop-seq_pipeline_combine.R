matrixc = readRDS(path.expand("~/GSE97930_CerebellarHem.RDS"))
matrixab = readRDS(path.expand("~/GSE97930_All.RDS"))

combinedmatrix = merge(matrixab, matrixc, by=0)

saveRDS(combinedmatrix, file = path.expand("~/GSE97930_All.RDS"))
