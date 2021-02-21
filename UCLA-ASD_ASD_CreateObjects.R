library(dplyr)
library(Seurat)
library(ggplot2)

##1BW
data_dir <- '/home/sheddn/UCLA-ASD/data/1BW_S1'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW1 = CreateSeuratObject(counts = expression_matrix)

##3BW
data_dir <- '/home/sheddn/UCLA-ASD/data/3BW_S3'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW3 = CreateSeuratObject(counts = expression_matrix)

##4BW
data_dir <- '/home/sheddn/UCLA-ASD/data/4BW_S4'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW3 = CreateSeuratObject(counts = expression_matrix)

##10BW
data_dir <- '/home/sheddn/UCLA-ASD/data/10BW_S13'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW10 = CreateSeuratObject(counts = expression_matrix)

##14BW
data_dir <- '/home/sheddn/UCLA-ASD/data/14BW_S17'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW14 = CreateSeuratObject(counts = expression_matrix)

##16BW
data_dir <- '/home/sheddn/UCLA-ASD/data/16BW_S2'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW16 = CreateSeuratObject(counts = expression_matrix)

##17BW
data_dir <- '/home/sheddn/UCLA-ASD/data/17BW_S3'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW17 = CreateSeuratObject(counts = expression_matrix)

##18BW
data_dir <- '/home/sheddn/UCLA-ASD/data/18BW_S4'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW18 = CreateSeuratObject(counts = expression_matrix)

##19BW
data_dir <- '/home/sheddn/UCLA-ASD/data/19BW_S5'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW19 = CreateSeuratObject(counts = expression_matrix)

##21BW
data_dir <- '/home/sheddn/UCLA-ASD/data/21BW_S7'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW21 = CreateSeuratObject(counts = expression_matrix)

##22BW
data_dir <- '/home/sheddn/UCLA-ASD/data/22BW_S8'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW22 = CreateSeuratObject(counts = expression_matrix)

##24BW
data_dir <- '/home/sheddn/UCLA-ASD/data/24BW_S10'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW24 = CreateSeuratObject(counts = expression_matrix)

##26BW
data_dir <- '/home/sheddn/UCLA-ASD/data/26BW_S2'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW26 = CreateSeuratObject(counts = expression_matrix)

##27W
data_dir <- '/home/sheddn/UCLA-ASD/data/27BW_S3'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW27 = CreateSeuratObject(counts = expression_matrix)

##28W
data_dir <- '/home/sheddn/UCLA-ASD/data/28BW_S4'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW28 = CreateSeuratObject(counts = expression_matrix)

##28W
data_dir <- '/home/sheddn/UCLA-ASD/data/29BW_S5'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW29 = CreateSeuratObject(counts = expression_matrix)

##31W
data_dir <- '/home/sheddn/UCLA-ASD/data/31BW_S1'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW31 = CreateSeuratObject(counts = expression_matrix)

##32W
data_dir <- '/home/sheddn/UCLA-ASD/data/32BW_S2'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW32 = CreateSeuratObject(counts = expression_matrix)

##34W
data_dir <- '/home/sheddn/UCLA-ASD/data/34BW_S4'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW34 = CreateSeuratObject(counts = expression_matrix)

##37W
data_dir <- '/home/sheddn/UCLA-ASD/data/37BW_S7'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW37 = CreateSeuratObject(counts = expression_matrix)

##38W
data_dir <- '/home/sheddn/UCLA-ASD/data/38BW_S8'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW38 = CreateSeuratObject(counts = expression_matrix)

##39W
data_dir <- '/home/sheddn/UCLA-ASD/data/39BW_S9'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW39 = CreateSeuratObject(counts = expression_matrix)

##42W
data_dir <- '/home/sheddn/UCLA-ASD/data/42BW_S12'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW42 = CreateSeuratObject(counts = expression_matrix)

##44W
data_dir <- '/home/sheddn/UCLA-ASD/data/44BW_S14'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW44 = CreateSeuratObject(counts = expression_matrix)

##45W
data_dir <- '/home/sheddn/UCLA-ASD/data/45BW_S15'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW45 = CreateSeuratObject(counts = expression_matrix)

##48W
data_dir <- '/home/sheddn/UCLA-ASD/data/48BW_S3'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW48 = CreateSeuratObject(counts = expression_matrix)

##49W
data_dir <- '/home/sheddn/UCLA-ASD/data/49BW_S4'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW49 = CreateSeuratObject(counts = expression_matrix)

##50W
data_dir <- '/home/sheddn/UCLA-ASD/data/50BW_S5'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW50 = CreateSeuratObject(counts = expression_matrix)

##54W
data_dir <- '/home/sheddn/UCLA-ASD/data/54BW_S11'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW54 = CreateSeuratObject(counts = expression_matrix)

##59W
data_dir <- '/home/sheddn/UCLA-ASD/data/59BW_S29'
list.files(data_dir) # Should show barcodes.tsv, genes.tsv, and matrix.mtx
expression_matrix <- Read10X(data.dir = data_dir)
BW59 = CreateSeuratObject(counts = expression_matrix)

## Missing 7, 9, 60, 65, 66, 67, 68, 70, 72, 73, 75, 76

ASD <- merge(BW1, y=c(BW3,BW4,BW10,BW14,BW16,BW17,BW18,BW19,BW21,BW22,BW24,BW26,BW27,BW28,
                      BW29,BW31,BW32,BW34,BW37,BW38,BW39,BW42,BW44,BW45,BW48,BW49,BW54,BW59),
             add.cell.ids=c('BW3','BW4','BW10','BW14','BW16','BW17','BW18','BW19','BW21','BW22','BW24','BW26','BW27','BW28',
                            'BW29','BW31','BW32','BW34','BW37','BW38','BW39','BW42','BW44','BW45','BW48','BW49','BW54','BW59'),
             , project = "UCLA-ASD")

saveRDS(ASD, '/home/sheddn/UCLA-ASD/data/ASD_SampleLabels')

ASD_forcombined <- merge(BW1, y=c(BW3,BW4,BW10,BW14,BW16,BW17,BW18,BW19,BW21,BW22,BW24,BW26,BW27,BW28,
                      BW29,BW31,BW32,BW34,BW37,BW38,BW39,BW42,BW44,BW45,BW48,BW49,BW54,BW59),
             add.cell.ids=c('ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD',
                            'ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD','ASD'),
             , project = "UCLA-ASD")

saveRDS(ASD, '/home/sheddn/UCLA-ASD/data/ASD_GroupLabel')

