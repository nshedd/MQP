matrix = read.table("~/subset.tsv", header=TRUE, row.names=1, sep='\t')

print("finsihed loading data")

saveRDS(matrix, file=path.expand("~/GRCh38-rDHSs_GW17_Cortex_50k_fragments_subset2.rds"))

print("finished saving RDS")
