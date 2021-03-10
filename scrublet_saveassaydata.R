ASD = readRDS('/data/rusers/sheddn/UCLA-ASD/data/ASD_BA4.6')

write.table(as.matrix(GetAssayData(object = seurat.obj, slot = "counts")), 
            "/data/rusers/sheddn/UCLA-ASD/data/CTL_BA4.6_counts.txt", 
            sep = '\t', row.names = T, col.names = T, quote = F)
