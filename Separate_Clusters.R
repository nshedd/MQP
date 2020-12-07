genes = read.table(path.expand("~/temporal_marker_genes.txt"), header=TRUE, row.names=1, sep="\t")

cluster1_genes = genes$name[genes$group == 1]
cluster2_genes = genes$name[genes$group == 2]
cluster3_genes = genes$name[genes$group == 3]
cluster4_genes = genes$name[genes$group == 4]
cluster5_genes = genes$name[genes$group == 5]
cluster6_genes = genes$name[genes$group == 6]
cluster7_genes = genes$name[genes$group == 7]
cluster8_genes = genes$name[genes$group == 8]
cluster9_genes = genes$name[genes$group == 9]
cluster10_genes = genes$name[genes$group == 10]
cluster11_genes = genes$name[genes$group == 11]

write.table(cluster1_genes, path.expand("~/temporal_genes/cluster1.txt"), sep = "\n")
write.table(cluster2_genes, path.expand("~/temporal_genes/cluster2.txt"), sep = "\n")
write.table(cluster3_genes, path.expand("~/temporal_genes/cluster3.txt"), sep = "\n")
write.table(cluster4_genes, path.expand("~/temporal_genes/cluster4.txt"), sep = "\n")
write.table(cluster5_genes, path.expand("~/temporal_genes/cluster5.txt"), sep = "\n")
write.table(cluster6_genes, path.expand("~/temporal_genes/cluster6.txt"), sep = "\n")
write.table(cluster7_genes, path.expand("~/temporal_genes/cluster7.txt"), sep = "\n")
write.table(cluster8_genes, path.expand("~/temporal_genes/cluster8.txt"), sep = "\n")
write.table(cluster9_genes, path.expand("~/temporal_genes/cluster9.txt"), sep = "\n")
write.table(cluster10_genes, path.expand("~/temporal_genes/cluster10.txt"), sep = "\n")
write.table(cluster11_genes, path.expand("~/temporal_genes/cluster11.txt"), sep = "\n")

