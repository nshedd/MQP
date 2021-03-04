library(dplyr)
library(Seurat)
library(ggplot2)

print('1BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/1BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW1 = CreateSeuratObject(counts = expression_matrix, project='1BW')

print('3BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/3BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW3 = CreateSeuratObject(counts = expression_matrix, project='3BW')

print('4BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/4BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW4 = CreateSeuratObject(counts = expression_matrix, project='4BW')

print('7BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/7BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW7 = CreateSeuratObject(counts = expression_matrix, project='7BW')

print('9BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/9BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW9 = CreateSeuratObject(counts = expression_matrix, project='9BW')

print('10BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/10BW_S13/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW10 = CreateSeuratObject(counts = expression_matrix, project='10BW')

print('14BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/14BW_S17/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW14 = CreateSeuratObject(counts = expression_matrix, project='14BW')

print('16BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/16BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW16 = CreateSeuratObject(counts = expression_matrix, project='16BW')

print('17BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/17BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW17 = CreateSeuratObject(counts = expression_matrix, project='17BW')

print('18BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/18BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW18 = CreateSeuratObject(counts = expression_matrix, project='18BW')

print('19BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/19BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW19 = CreateSeuratObject(counts = expression_matrix, project='19BW')

print('21BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/21BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW21 = CreateSeuratObject(counts = expression_matrix, project='21BW')

print('22BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/22BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW22 = CreateSeuratObject(counts = expression_matrix, project='22BW')

print('24BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/24BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW24 = CreateSeuratObject(counts = expression_matrix, project='24BW')

print('26BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/26BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW26 = CreateSeuratObject(counts = expression_matrix, project='26BW')

print('27BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/27BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW27 = CreateSeuratObject(counts = expression_matrix, project='27BW')

print('28BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/28BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW28 = CreateSeuratObject(counts = expression_matrix, project='28BW')

print('29BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/29BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW29 = CreateSeuratObject(counts = expression_matrix, project='29BW')

print('31BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/31BW_S1/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW31 = CreateSeuratObject(counts = expression_matrix, project='31BW')

print('32BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/32BW_S2/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW32 = CreateSeuratObject(counts = expression_matrix, project='32BW')

print('34BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/34BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW34 = CreateSeuratObject(counts = expression_matrix, project='34BW')

print('37BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/37BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW37 = CreateSeuratObject(counts = expression_matrix, project='37BW')

print('38BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/38BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW38 = CreateSeuratObject(counts = expression_matrix, project='38BW')

print('39BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/39BW_S9/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW39 = CreateSeuratObject(counts = expression_matrix, project='39BW')

print('42BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/42BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW42 = CreateSeuratObject(counts = expression_matrix, project='42BW')

print('44BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/44BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW44 = CreateSeuratObject(counts = expression_matrix, project='44BW')

print('45BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/45BW_S15/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW45 = CreateSeuratObject(counts = expression_matrix, project='45BW')

print('48BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/48BW_S3/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW48 = CreateSeuratObject(counts = expression_matrix, project='48BW')

print('49BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/49BW_S4/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW49 = CreateSeuratObject(counts = expression_matrix, project='49BW')

print('50BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/50BW_S5/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW50 = CreateSeuratObject(counts = expression_matrix, project='50BW')

print('54BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/54BW_S11/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW54 = CreateSeuratObject(counts = expression_matrix, project='54BW')

print('59BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/59BW_S29/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW59 = CreateSeuratObject(counts = expression_matrix, project='59BW')

print('60BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/60BW_S30/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW60 = CreateSeuratObject(counts = expression_matrix, project='60BW')

print('65BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/65BW_S7/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW65 = CreateSeuratObject(counts = expression_matrix, project='65BW')

print('66BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/66BW_S8/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW66 = CreateSeuratObject(counts = expression_matrix, project='66BW')

print('67BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/67BW_S9/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW67 = CreateSeuratObject(counts = expression_matrix, project='67BW')

print('68BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/68BW_S10/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW68 = CreateSeuratObject(counts = expression_matrix, project='68BW')

print('70BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/70BW_S12/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW70 = CreateSeuratObject(counts = expression_matrix, project='70BW')

print('72BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/72BW_S14/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW72 = CreateSeuratObject(counts = expression_matrix, project='72BW')

print('73BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/73BW_S15/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW73 = CreateSeuratObject(counts = expression_matrix, project='73BW')

print('75BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/75BW_S17/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW75 = CreateSeuratObject(counts = expression_matrix, project='75BW')

print('76BW')
data_dir <- '/data/rusers/sheddn/UCLA-ASD/data/76BW_S18/Solo.out/GeneFull/filtered'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW76 = CreateSeuratObject(counts = expression_matrix, project='76BW')

ASD_BA4.6 <- merge(BW1, y=c(BW10,BW14,BW16,BW17,BW22,BW24,BW26,BW27,BW3,BW31,BW32,BW4,BW59,BW60,BW7,BW9),
                    add.cell.ids=c('BW1','BW10','BW14','BW16','BW17','BW22','BW24','BW26','BW27','BW3','BW31','BW32','BW4','BW59','BW60','BW7','BW9'),
                   project = "UCLA-ASD")

saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_BA4.6')

ASD_BA9 <- merge(BW18, y=c(BW19,BW21,BW28,BW29,BW34,BW37,BW38,BW39,BW42,BW44,BW45,BW48,BW49,BW50,BW54,BW65,BW66,BW67,BW68,BW70,BW72,BW73,BW75,BW76),
                    add.cell.ids=c('BW18','BW19','BW21','BW28','BW29','BW34','BW37','BW38','BW39','BW42','BW44','BW45','BW48','BW49','BW50','BW54','BW65','BW66','BW67','BW68','BW70','BW72','BW73','BW75','BW76'),
                   project = "UCLA-ASD")

saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_BA9')


# ASD <- merge(BW1, y=c(BW3,BW4,BW7,BW9,BW10,BW14,BW16,BW17,BW18,BW19,BW21,BW22,BW24,BW26,BW27,BW28,BW29,BW31,BW32,BW34,
#                       BW37,BW38,BW39,BW42,BW44,BW45,BW48,BW49,BW54,BW59,BW60,BW65,BW66,BW67,BW68,BW70,BW72,BW73,BW75,BW76),
#              add.cell.ids=c('BW1','BW3','BW4','BW6','BW9','BW10','BW14','BW16','BW17','BW18','BW19','BW21','BW22','BW24',
#                             'BW26','BW27','BW28','BW29','BW31','BW32','BW34','BW37','BW38','BW39','BW42','BW44','BW45',
#                             'BW48','BW49','BW54','BW59','BW60','BW65','BW66','BW67','BW68','BW70','BW72','BW73','BW75','BW76'),
#              project = "UCLA-ASD")
# 
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_SampleLabels')
# 
# new.cluster.ids <- c('ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD',
#                      'ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD',
#                      'ASD','ASD','ASD','ASD','ASD','ASD','ASD')
# names(new.cluster.ids) <- levels(ASD)
# ASD <- RenameIdents(ASD, new.cluster.ids)
# 
# saveRDS(ASD, '/data/rusers/sheddn/UCLA-ASD/data/ASD_GroupLabel')

