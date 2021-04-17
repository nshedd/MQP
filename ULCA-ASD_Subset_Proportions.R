library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA4.6_WithDEGs.RDS')

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

write.table(BA4.6_prop_df, '/data/rusers/sheddn/UCLA-ASD/subset/data/CellTypeProportions_BA4.6.txt', sep=',')




BA9 <- readRDS('/data/rusers/sheddn/UCLA-ASD/subset/data/combined_BA9_WithDEGs.RDS')

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

write.table(BA9_prop_df, '/data/rusers/sheddn/UCLA-ASD/subset/data/CellTypeProportions_BA9.txt', sep=',')
