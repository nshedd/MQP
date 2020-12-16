install.packages("sqldf")


path1a = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")
path1b = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")
path1c = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrixa = read.table(path1a, header=TRUE, row.names=1)
matrixb = read.table(path1b, header=TRUE, row.names=1)
matrixc = read.table(path1c, header=TRUE, row.names=1)

saveRDS(matrixa, path.expand("~/GSE97930_VisualCortex.RDS"))
saveRDS(matrixb, path.expand("~/GSE97930_FrontalCortex.RDS"))
saveRDS(matrixc, path.expand("~/GSE97930_CerebellarHem.RDS"))

library(sqldf)

matrixa = file(matrixa)
matrixb = file(matrixb)
matrixc = file(matrixc)

matrix = sqldf("select * from matrixa union select * from matrixb union select * from matrixc", dbname = tempfile())

saveRDS(matrix, file = path.expand("~/GSE97930_All_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"))
