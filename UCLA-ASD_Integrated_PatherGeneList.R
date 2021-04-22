Ast <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Ast.txt")
Ast_genes <- row.names(Ast)
write.table(Ast_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Ast.txt", sep=',')

Ex <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Ex.txt")
Ex_genes <- row.names(Ex)
write.table(Ex_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Ex.txt", sep=',')

In <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_In.txt")
In_genes <- row.names(In)
write.table(In_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_In.txt", sep=',')

Mic <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Mic.txt")
Mic_genes <- row.names(Mic)
write.table(Mic_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Mic.txt", sep=',')

Oli <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_Oli.txt")
Oli_genes <- row.names(Oli)
write.table(Oli_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_Oli.txt", sep=',')

OPC <- read.table("/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-FDR_OPC.txt")
OPC_genes <- row.names(OPC)
write.table(OPC_genes, "/data/rusers/sheddn/UCLA-ASD/plots/DEG-nomito/BA4.6_DEGs-genelist_OPC.txt", sep=',')
