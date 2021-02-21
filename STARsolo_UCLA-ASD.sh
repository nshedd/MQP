work_path1=/home/sheddn/UCLA-ASD/data
work_path2=/data/projects/psychencode/TR-SingleCell/UCLA-ASD/PEC_syn18898607/scRNAseq

  
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/SCZ-BP/ \
--readFilesIn <(gunzip -c ${work_path2}/10BW_S13_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/10BW_S13_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/10BW_S13_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/10BW_S13_L002_R1_001.fastq.gz) \
--soloType Droplet \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 16 \
--soloBarcodeReadLength 0 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures Gene \
--outFileNamePrefix ${work_path1}/10BW_S13/
