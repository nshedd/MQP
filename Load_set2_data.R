matrix = read.table("/data/zusers/pratth/sc/atac/GSM3722075_PBMC_Rep3_fragments.tsv.gz.rDHS.matrix.tsv", header=TRUE, sep="/t", row.names=1)

print("finsihed loading data")

saveRDS(matrix, file=path.expand("~/GSM3722075_PBMC_Rep3_fragments.rds"))

print("finished saving RDS")