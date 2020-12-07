genes = read.table(path.expand("~/GSE97930_visualcortex_differentiallyexpressed.txt"), header=TRUE, row.names=1, sep="\t")

cluster1_genes = genes$gene[genes$cluster == '0']
cluster1_genes = genes$gene[genes$cluster == '1']
cluster2_genes = genes$gene[genes$cluster == '2']
cluster3_genes = genes$gene[genes$cluster == '3']
cluster4_genes = genes$gene[genes$cluster == '4']
cluster5_genes = genes$gene[genes$cluster == '5']
cluster6_genes = genes$gene[genes$cluster == '6']
cluster7_genes = genes$gene[genes$cluster == '7']
cluster8_genes = genes$gene[genes$cluster == '8']
cluster9_genes = genes$gene[genes$cluster == '9']
cluster10_genes = genes$gene[genes$cluster == '10']
cluster11_genes = genes$gene[genes$cluster == '11']


write.table(cluster0_genes, path.expand("~/GSE97930_visualcortex_genes/cluster0.txt"), sep = "\n", row.names=F)
write.table(cluster1_genes, path.expand("~/GSE97930_visualcortex_genes/cluster1.txt"), sep = "\n", row.names=F)
write.table(cluster2_genes, path.expand("~/GSE97930_visualcortex_genes/cluster2.txt"), sep = "\n", row.names=F)
write.table(cluster3_genes, path.expand("~/GSE97930_visualcortex_genes/cluster3.txt"), sep = "\n", row.names=F)
write.table(cluster4_genes, path.expand("~/GSE97930_visualcortex_genes/cluster4.txt"), sep = "\n", row.names=F)
write.table(cluster5_genes, path.expand("~/GSE97930_visualcortex_genes/cluster5.txt"), sep = "\n", row.names=F)
write.table(cluster6_genes, path.expand("~/GSE97930_visualcortex_genes/cluster6.txt"), sep = "\n", row.names=F)
write.table(cluster7_genes, path.expand("~/GSE97930_visualcortex_genes/cluster7.txt"), sep = "\n", row.names=F)
write.table(cluster8_genes, path.expand("~/GSE97930_visualcortex_genes/cluster8.txt"), sep = "\n", row.names=F)
write.table(cluster9_genes, path.expand("~/GSE97930_visualcortex_genes/cluster9.txt"), sep = "\n", row.names=F)
write.table(cluster10_genes, path.expand("~/GSE97930_visualcortex_genes/cluster10.txt"), sep = "\n", row.names=F)
write.table(cluster11_genes, path.expand("~/GSE97930_visualcortex_genes/cluster11.txt"), sep = "\n", row.names=F)

