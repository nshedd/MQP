library(dplyr)
library(Seurat)
library(ggplot2)

print('2BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/2BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW2 = CreateSeuratObject(counts = expression_matrix, project='2BW')

print('5BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/5BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW5 = CreateSeuratObject(counts = expression_matrix, project='5BW')

print('6BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/6BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW6 = CreateSeuratObject(counts = expression_matrix, project='6BW')

print('8BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/8BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW8 = CreateSeuratObject(counts = expression_matrix, project='8BW')

print('11BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/11BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW11 = CreateSeuratObject(counts = expression_matrix, project='11BW')

print('12BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/12BW_S15/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW12 = CreateSeuratObject(counts = expression_matrix, project='12BW')

print('13BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/13BW_S16/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW13 = CreateSeuratObject(counts = expression_matrix, project='13BW')

print('15BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/15BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW15 = CreateSeuratObject(counts = expression_matrix, project='15BW')

print('20BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/20BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW20 = CreateSeuratObject(counts = expression_matrix, project='20BW')

print('23BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/23BW_S9/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW23 = CreateSeuratObject(counts = expression_matrix, project='23BW')

print('25BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/25BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW25 = CreateSeuratObject(counts = expression_matrix, project='25BW')

print('30BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/30BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW30 = CreateSeuratObject(counts = expression_matrix, project='30BW')

print('33BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/33BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW33 = CreateSeuratObject(counts = expression_matrix, project='33BW')

print('35BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/35BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW35 = CreateSeuratObject(counts = expression_matrix, project='35BW')

print('40BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/40BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW40 = CreateSeuratObject(counts = expression_matrix, project='40BW')

print('41BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/41BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW41 = CreateSeuratObject(counts = expression_matrix, project='41BW')

print('43BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/43BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW43 = CreateSeuratObject(counts = expression_matrix, project='43BW')

print('46BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/46BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW46 = CreateSeuratObject(counts = expression_matrix, project='46BW')

print('47BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/47BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW47 = CreateSeuratObject(counts = expression_matrix, project='47BW')

print('51BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/51BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW51 = CreateSeuratObject(counts = expression_matrix, project='51BW')

print('52BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/52BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW52 = CreateSeuratObject(counts = expression_matrix, project='52BW')

print('53BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/53BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW53 = CreateSeuratObject(counts = expression_matrix, project='53BW')

print('55BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/55BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW55 = CreateSeuratObject(counts = expression_matrix, project='55BW')

print('56BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/56BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW56 = CreateSeuratObject(counts = expression_matrix, project='56BW')

print('57BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/57BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW57 = CreateSeuratObject(counts = expression_matrix, project='57BW')

print('58BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/58BW_S28/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW58 = CreateSeuratObject(counts = expression_matrix, project='58BW')

print('61BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/61BW_S31/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW61 = CreateSeuratObject(counts = expression_matrix, project='61BW')

print('62BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/62BW_S32/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW62 = CreateSeuratObject(counts = expression_matrix, project='62BW')

print('63BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/63BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW63 = CreateSeuratObject(counts = expression_matrix, project='63BW')

print('64BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/64BW_S6/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW64 = CreateSeuratObject(counts = expression_matrix, project='64BW')

print('69BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/69BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW69 = CreateSeuratObject(counts = expression_matrix, project='69BW')

print('71BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/71BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW71 = CreateSeuratObject(counts = expression_matrix, project='71BW')

print('74BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/74BW_S16/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW74 = CreateSeuratObject(counts = expression_matrix, project='74BW')

print('77BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/77BW_S19/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW77 = CreateSeuratObject(counts = expression_matrix, project='77BW')

print('78BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/78BW_S20/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW78 = CreateSeuratObject(counts = expression_matrix, project='78BW')

# CTL_BA4.6 <- merge(BW11, y=c(BW12,BW13,BW15,BW2,BW23,BW25,BW30,BW35,BW46,BW5,BW6,BW8),
#                     add.cell.ids=c('BW11','BW12','BW13','BW15','BW2','BW23','BW25','BW30','BW35','BW46','BW5','BW6','BW8'),
#                    project = "UCLA-ASD")
# 
# saveRDS(CTL_BA4.6, '/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6')

CTL_BA9 <- merge(BW20, y=c(BW33,BW40,BW41,BW43,BW47,BW51,BW52,BW53,BW55,BW56,BW57,BW61,BW62,BW63,BW64,BW69,BW71,BW74,BW77,BW78),
                    add.cell.ids=c('BW20','BW33','BW40','BW41','BW43','BW47','BW51','BW52','BW53','BW55','BW56','BW57','BW61','BW62','BW63','BW64','BW69','BW71','BW74','BW77','BW78'),
                   project = "UCLA-ASD")

saveRDS(CTL_BA9, '/data/rusers/sheddn/UCLA-ASD/data/CTL_BA9')


# CTL <- merge(BW2, y=c(BW5,BW6,BW8,BW11,BW12,BW13,BW15,BW20,BW23,BW25,BW30,BW33,BW35,BW40,BW41,BW43,BW46,BW47,
#                       BW51,BW52,BW53,BW55,BW56,BW57,BW58,BW61,BW62,BW63,BW64,BW69,BW71,BW74,BW77,BW78),
#              add.cell.ids=c('BW2','BW5','BW6','BW8','BW11','BW12','BW13','BW15','BW20','BW23','BW25','BW30','BW33','BW35','BW40','BW41','BW43','BW46',
#                             'BW47','BW51','BW52','BW53','BW55','BW56','BW57','BW58','BW61','BW62','BW63','BW64','BW69','BW71','BW74','BW77','BW78'),
#              project = "UCLA-ASD")
# 
# 
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_SampleLabels')
# 
# new.cluster.ids <- c('CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL',
#                      'CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL','CTL', 'CTL')
# names(new.cluster.ids) <- levels(CTL)
# CTL <- RenameIdents(CTL, new.cluster.ids)
# 
# saveRDS(CTL, '/data/rusers/sheddn/UCLA-ASD/data/CTL_GroupLabel')


