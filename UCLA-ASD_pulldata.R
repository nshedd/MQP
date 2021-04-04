library(Seurat)

BA4.6 <- readRDS('/data/rusers/sheddn/UCLA-ASD/data/combined_BA4.6_WithDEGs.RDS')

data.expression = AverageExpression(data, slot="counts")
data.expresssion = as.matrix(data.expression[["RNA"]])
head(data.expression)
q()


group = c('CTL','ASD')
sample = unique(BA4.6$Sample)
cluster = unique(Idents(BA4.6))

expression.0 <- expression.clusters
expression.1 <- expression.clusters
expression.2 <- expression.clusters
expression.3 <- expression.clusters
expression.4 <- expression.clusters           
expression.5 <- expression.clusters               
expression.6 <- expression.clusters
expression.7 <- expression.clusters
expression.8 <- expression.clusters
expression.9 <- expression.clusters
expression.10 <- expression.clusters
expression.11 <- expression.clusters
expression.12 <- expression.clusters
expression.13 <- expression.clusters
expression.14 <- expression.clusters              
expression.15 <- expression.clusters               
expression.16 <- expression.clusters
expression.17 <- expression.clusters
expression.18 <- expression.clusters
expression.19 <- expression.clusters              
expression.20 <- expression.clusters   
expression.21 <- expression.clusters
expression.22 <- expression.clusters
expression.23 <- expression.clusters
expression.24 <- expression.clusters              
expression.25 <- expression.clusters               
expression.26 <- expression.clusters
expression.27 <- expression.clusters
expression.28 <- expression.clusters
expression.29 <- expression.clusters              
expression.30 <- expression.clusters
expression.31 <- expression.clusters
expression.32 <- expression.clusters
expression.33 <- expression.clusters
expression.34 <- expression.clusters              
expression.35 <- expression.clusters               
expression.36 <- expression.clusters
expression.37 <- expression.clusters
expression.38 <- expression.clusters
              
               
for (s in sample) {
  data <- subset(BA4.6, subset = Sample==s)
  data.expression = AverageExpression(data, slot="counts")
  data.expresssion = as.matrix(data.expression[["RNA"]])
  expression.0 <- merge(expression.0, data.expression$C0, all=TRUE)
  expression.1 <- merge(expression.1, data.expression$C1, all=TRUE)
  expression.2 <- merge(expression.2, data.expression$C2, all=TRUE)
  expression.3 <- merge(expression.3, data.expression$C3, all=TRUE)
  expression.4 <- merge(expression.4, data.expression$C4, all=TRUE)           
  expression.5 <- merge(expression.5, data.expression$C5, all=TRUE)              
  expression.6 <- merge(expression.6, data.expression$C6, all=TRUE)
  expression.7 <- merge(expression.7, data.expression$C7, all=TRUE)
  expression.8 <- merge(expression.8, data.expression$C8, all=TRUE)
  expression.9 <- merge(expression.9, data.expression$C9, all=TRUE)
  expression.10 <- merge(expression.10, data.expression$C10, all=TRUE)
  expression.11 <- merge(expression.11, data.expression$C11, all=TRUE)
  expression.12 <- merge(expression.12, data.expression$C12, all=TRUE)
  expression.13 <- merge(expression.13, data.expression$C13, all=TRUE)
  expression.14 <- merge(expression.14, data.expression$C14, all=TRUE)              
  expression.15 <- merge(expression.15, data.expression$C15, all=TRUE)               
  expression.16 <- merge(expression.16, data.expression$C16, all=TRUE)
  expression.17 <- merge(expression.17, data.expression$C17, all=TRUE)
  expression.18 <- merge(expression.18, data.expression$C18, all=TRUE)
  expression.19 <- merge(expression.19, data.expression$C19, all=TRUE)             
  expression.20 <- merge(expression.20, data.expression$C20, all=TRUE)   
  expression.21 <- merge(expression.21, data.expression$C21, all=TRUE)
  expression.22 <- merge(expression.22, data.expression$C22, all=TRUE)
  expression.23 <- merge(expression.23, data.expression$C23, all=TRUE)
  expression.24 <- merge(expression.24, data.expression$C24, all=TRUE)           
  expression.25 <- merge(expression.25, data.expression$C25, all=TRUE)            
  expression.26 <- merge(expression.26, data.expression$C26, all=TRUE)
  expression.27 <- merge(expression.27, data.expression$C27, all=TRUE)
  expression.28 <- merge(expression.28, data.expression$C28, all=TRUE)
  expression.29 <- merge(expression.29, data.expression$C29, all=TRUE)  
  expression.30 <- merge(expression.30, data.expression$C30, all=TRUE)
  expression.31 <- merge(expression.31, data.expression$C31, all=TRUE)
  expression.32 <- merge(expression.32, data.expression$C32, all=TRUE)
  expression.33 <- merge(expression.33, data.expression$C33, all=TRUE)
  expression.34 <- merge(expression.34, data.expression$C34, all=TRUE)        
  expression.35 <- merge(expression.35, data.expression$C35, all=TRUE)             
  expression.36 <- merge(expression.36, data.expression$C36, all=TRUE)
  expression.37 <- merge(expression.37, data.expression$C37, all=TRUE)
  expression.38 <- merge(expression.38, data.expression$C38, all=TRUE)
}
              
head(expression.0)
