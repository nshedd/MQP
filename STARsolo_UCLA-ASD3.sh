work_path1=/data/rusers/sheddn/UCLA-ASD/data
work_path2=/data/projects/psychencode/TR-SingleCell/UCLA-ASD/PEC_syn18898607/scRNAseq

labels3="31BW_S1 32B_S2 33BW_S3 34BW_S4 35BW_S5 36BW_S6 37BW_S7 38BW_S8 39BW_S9 40BW_S10 41BW_S11 42BW_S12 43BW_S13 44BW_S14 45BW_S15"

for f in $labels3
do
	echo ${f}
	--genomeDir ~/UCLA-ASD/ \
	--readFilesIn <(gunzip -c ${work_path2}/${f}_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/${f}_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L004_R1_001.fastq.gz), <(gunzip -c ${work_path2}/${f}_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/${f}_L006_R1_001.fastq.gz) \
	--soloType CB_UMI_Simple \
	--soloCBwhitelist ~/3M-february-2018.txt \
	--runThreadN 12 \
	--soloBarcodeReadLength 1 \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
	--soloFeatures GeneFull \
	--outFileNamePrefix ${work_path1}/${f}/
done
