install.packages("sqldf")
library(sqldf)

path1a = path.expand("~/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")
path1b = path.expand("~/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz")
path1c = path.expand("~/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt")

matrixa = file(path1a)
matrixb = file(path1b)
matrixc = file(path1c)

matrix = sqldf("select * from matrixa union select * from matrixb union select * from matrixc", dbname = tempfile())

saveRDS(matrix, file = path.expand("~/GSE97930_All_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt.gz"))