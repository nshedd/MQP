work_path1=/data/zusers/sheddn
work_path2=/home/sheddn/UCLA-ASD/PEC_syn18898607/scRNAseq

  
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/SCZ-BP/ \
--readFilesIn ${work_path2}/12BW_S15_L001_R1_001.fastq ${work_path2}/12BW_S15_L001_R2_001.fastq \
--soloType Droplet \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 16 \
--soloBarcodeReadLength 0 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures Gene \
--outFileNamePrefix ${work_path2}/12BW_S15/