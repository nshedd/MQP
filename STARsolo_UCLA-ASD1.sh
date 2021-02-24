work_path1=/data/rusers/sheddn/UCLA-ASD/data
work_path2=/data/projects/psychencode/TR-SingleCell/UCLA-ASD/PEC_syn18898607/scRNAseq

labels="1BW_S1 2BW_S2 3BW_S3 4BW_S4 5BW_S5 6BW_S6 7BW_S10 8BW_S11 9BW_S12 10BW_S13 11BW_S14 12BW_S15 13BW_S16 14BW_S17"

for f in $labels
do
	echo ${f}
	# Run STARsolo 2.7.7a
	~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
	--genomeDir ~/UCLA-ASD/ \
	--readFilesIn <(gunzip -c ${work_path2}/${f}_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/${f}_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L002_R1_001.fastq.gz) \
	--soloType CB_UMI_Simple \
	--soloCBwhitelist ~/737K-august-2016.txt \
	--runThreadN 12 \
	--soloBarcodeReadLength 1 \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
	--soloFeatures GeneFull \
	--outFileNamePrefix ${work_path1}/${f}/
done
