matrix = read.table("/data/zusers/pratth/sc/atac/PEC/GRCh38-rDHSs_GW17_Cortex.aggregate.50k.tsv", header=TRUE, row.names=1, sep='\t')

print("finsihed loading data")

saveRDS(matrix, file=path.expand("~/GRCh38-rDHSs_GW17_Cortex_50k_fragments.rds"))

print("finished saving RDS")
