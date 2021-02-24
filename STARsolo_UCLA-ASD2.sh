work_path1=/data/rusers/sheddn/UCLA-ASD/data
work_path2=/data/projects/psychencode/TR-SingleCell/UCLA-ASD/PEC_syn18898607/scRNAseq

labels2="15BW_S1 16BW_S2 17BW_S3 18BW_S4 19BW_S5 20BW_S6 21BW_S7 22BW_S8 23BW_S9 24BW_S10 25BW_S1 26BW_S2 27BW_S3 28BW_S4 29BW_S5 30BW_S6"

for f in $labels2
do
	echo ${f}
	# Run STARsolo 2.7.7a
	~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
	--genomeDir ~/UCLA-ASD/ \
	--readFilesIn <(gunzip -c ${work_path2}/${f}_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/${f}_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L002_R1_001.fastq.gz) \
	--soloType CBUMI_Simple \
	--soloCBwhitelist ~/3M-february-2018.txt \
	--runThreadN 12 \
	--soloBarcodeReadLength 1 \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
	--soloFeatures GeneFull \
	--outFileNamePrefix ${work_path1}/${f}/
done
