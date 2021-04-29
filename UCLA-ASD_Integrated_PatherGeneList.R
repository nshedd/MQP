Ast <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Ast.txt")
Ast_genes <- row.names(Ast)
print(Ast_genes)
write.table(Ast_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Ast.txt")

Ex <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Ex.txt")
Ex_genes <- row.names(Ex)
print(Ex_genes)
write.table(Ex_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Ex.txt")

In <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_In.txt")
In_genes <- row.names(In)
print(In_genes)
write.table(In_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_In.txt")

Mic <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Mic.txt")
Mic_genes <- row.names(Mic)
print(Mic_genes)
write.table(Mic_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Mic.txt")

Oli <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Oli.txt")
Oli_genes <- row.names(Oli)
print(Oli_genes)
write.table(Oli_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Oli.txt")

OPC <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_OPC.txt")
OPC_genes <- row.names(OPC)
print(OPC_genes)
write.table(OPC_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_OPC.txt")
