library(dplyr)
library(Seurat)
library(ggplot2)


## Read and relabel clusters and plot
BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA9_WithDEGs.RDS')

p1 <- DimPlot(BA9, reduction = "umap", group.by = "Group")
p2 <- DimPlot(BA9, reduction = "umap", label = TRUE)
p3 <- p1 + p2
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_combined.png', width = 16, height = 7)

new_ids = c("Ex1","Ex2","Ex3","In1","Ast1","OPC1","Mic","Ex4","In2","Ex5",
            "In3","Ex6","Ex7","In4","Ex8","Oli","Ast2","Ex9","Ex10","Ex11",
            "Ex11","Ex12","Ex13","Ex14","Ast3","In5","In6","Ex15","Ex16","In7",
            "OPC2","Ex17")

new.cluster.ids <- new_ids
names(new.cluster.ids) <- levels(BA9)
BA9 <- RenameIdents(BA9, new.cluster.ids)

DimPlot(BA9, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel_bycluster_relabeled.png', width = 8, height = 7)

DimPlot(BA9, reduction = "umap", label=TRUE, split.by = "Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA9_integrated_SingleRlabel_bycluster_separate_relabeled.png', width = 16, height = 7)

## By cluster proportion analysis
ASD_BA9 <- subset(BA9, subset = Group == "ASD")

ASD_BA9_prop_manual_list = Idents(ASD_BA9)
ASD_BA9_prop_manual_table = table(ASD_BA9_prop_manual_list)

ASD_BA9_prop_df = ASD_BA9_prop_manual_table %>% as.data.frame
colnames(ASD_BA9_prop_df) <- c('Cell_Type', 'FreqASD')

ASD_BA9_sum = sum(ASD_BA9_prop_df$FreqASD)
ASD_BA9_prop_df$FreqASD = ASD_BA9_prop_df$FreqASD/ASD_BA9_sum

for (i in unique(ASD_BA9$Sample)) {
  data <- subset(ASD_BA9, subset = Sample == i)
  
  lab = paste("ASD", i)
  
  ASD_BA9_prop_manual_list = Idents(data)
  ASD_BA9_prop_manual_table = table(ASD_BA9_prop_manual_list)
  
  df = ASD_BA9_prop_manual_table %>% as.data.frame
  colnames(df) <- c('Cell_Type', lab)
  
  ASD_BA9_sum = sum(df[lab])
  df[lab] = df[lab]/ASD_BA9_sum
  
  ASD_BA9_prop_df = merge(x = ASD_BA9_prop_df, y = df, by = 'Cell_Type', all = TRUE)
}


CTL_BA9 <- subset(BA9, subset = Group == "CTL")

CTL_BA9_prop_manual_list = Idents(CTL_BA9)
CTL_BA9_prop_manual_table = table(CTL_BA9_prop_manual_list)

CTL_BA9_prop_df = CTL_BA9_prop_manual_table %>% as.data.frame
colnames(CTL_BA9_prop_df) <- c('Cell_Type', 'FreqCTL')

CTL_BA9_sum = sum(CTL_BA9_prop_df$FreqCTL)
CTL_BA9_prop_df$FreqCTL = CTL_BA9_prop_df$FreqCTL/CTL_BA9_sum

for (i in unique(CTL_BA9$Sample)) {
  data <- subset(CTL_BA9, subset = Sample == i)
  
  lab = paste("CTL", i)
  
  CTL_BA9_prop_manual_list = Idents(data)
  CTL_BA9_prop_manual_table = table(CTL_BA9_prop_manual_list)
  
  df = CTL_BA9_prop_manual_table %>% as.data.frame
  colnames(df) <- c('Cell_Type', lab)
  
  CTL_BA9_sum = sum(df[lab])
  df[lab] = df[lab]/CTL_BA9_sum
  
  CTL_BA9_prop_df = merge(x = CTL_BA9_prop_df, y = df, by = 'Cell_Type', all = TRUE)
}

BA9_prop_df = merge(x = ASD_BA9_prop_df, y = CTL_BA9_prop_df, by = 'Cell_Type', all = TRUE)
print(BA9_prop_df)

write.table(BA9_prop_df, '/data/rusers/sheddn/UCLA-ASD/data/CellTypeProportions_SingleRbyclust_BA9.txt', sep=',')

q()


BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

new_ids = c("Ex1","Ex2","Ex3","In1","Ast1","Oli1","In2","OPC1","In3","Ex4",
            "Ex5","Ex6","Ex7","Ex8","Ex9","In4","Ex10","Ex11","In5","Ex12",
            "Ex13","Oli2","Ex14","Ast2","Ast3","In6","Ast4","In7","Mic1","In8",
            "Oli3","In9","Ex15","In10","In11","OPC2","In12","Ex16","Ast5","In13",
            "Ast6","OPC3","Mic2")


new.cluster.ids <- new_ids
names(new.cluster.ids) <- levels(BA4.6)
BA4.6 <- RenameIdents(BA4.6, new.cluster.ids)

DimPlot(BA4.6, label=TRUE, pt.size=0.5)
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel_bycluster_relabeled.png', width = 8, height = 7)

DimPlot(BA4.6, reduction = "umap", label=TRUE, split.by = "Group")
ggsave('/data/rusers/sheddn/UCLA-ASD/plots/UMAP_Harmony_BA4.6_integrated_SingleRlabel_bycluster_separate_relabeled.png', width = 16, height = 7)

## By cluster proportion analysis
ASD_BA4.6 <- subset(BA4.6, subset = Group == "ASD")

ASD_BA4.6_prop_manual_list = Idents(ASD_BA4.6)
ASD_BA4.6_prop_manual_table = table(ASD_BA4.6_prop_manual_list)

ASD_BA4.6_prop_df = ASD_BA4.6_prop_manual_table %>% as.data.frame
colnames(ASD_BA4.6_prop_df) <- c('Cell_Type', 'FreqASD')

ASD_BA4.6_sum = sum(ASD_BA4.6_prop_df$FreqASD)
ASD_BA4.6_prop_df$FreqASD = ASD_BA4.6_prop_df$FreqASD/ASD_BA4.6_sum

for (i in unique(ASD_BA4.6$Sample)) {
  data <- subset(ASD_BA4.6, subset = Sample == i)

  lab = paste("ASD", i)

  ASD_BA4.6_prop_manual_list = Idents(data)
  ASD_BA4.6_prop_manual_table = table(ASD_BA4.6_prop_manual_list)

  df = ASD_BA4.6_prop_manual_table %>% as.data.frame
  colnames(df) <- c('Cell_Type', lab)

  ASD_BA4.6_sum = sum(df[lab])
  df[lab] = df[lab]/ASD_BA4.6_sum

  ASD_BA4.6_prop_df = merge(x = ASD_BA4.6_prop_df, y = df, by = 'Cell_Type', all = TRUE)

}


CTL_BA4.6 <- subset(BA4.6, subset = Group == "CTL")

CTL_BA4.6_prop_manual_list = Idents(CTL_BA4.6)
CTL_BA4.6_prop_manual_table = table(CTL_BA4.6_prop_manual_list)

CTL_BA4.6_prop_df = CTL_BA4.6_prop_manual_table %>% as.data.frame
colnames(CTL_BA4.6_prop_df) <- c('Cell_Type', 'FreqCTL')

CTL_BA4.6_sum = sum(CTL_BA4.6_prop_df$FreqCTL)
CTL_BA4.6_prop_df$FreqCTL = CTL_BA4.6_prop_df$FreqCTL/CTL_BA4.6_sum

for (i in unique(CTL_BA4.6$Sample)) {
  data <- subset(CTL_BA4.6, subset = Sample == i)

  lab = paste("CTL", i)

  CTL_BA4.6_prop_manual_list = Idents(data)
  CTL_BA4.6_prop_manual_table = table(CTL_BA4.6_prop_manual_list)

  df = CTL_BA4.6_prop_manual_table %>% as.data.frame
  colnames(df) <- c('Cell_Type', lab)

  CTL_BA4.6_sum = sum(df[lab])
  df[lab] = df[lab]/CTL_BA4.6_sum

  CTL_BA4.6_prop_df = merge(x = CTL_BA4.6_prop_df, y = df, by = 'Cell_Type', all = TRUE)
}

BA4.6_prop_df = merge(x = ASD_BA4.6_prop_df, y = CTL_BA4.6_prop_df, by = 'Cell_Type', all = TRUE)
print(BA4.6_prop_df)

write.table(BA4.6_prop_df, '/data/rusers/sheddn/UCLA-ASD/data/CellTypeProportions_SingleRbyclust_BA4.6.txt', sep=',')

q()

