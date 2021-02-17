
work_path1=/data/zusers/sheddn
work_path2=/data/zusers/pratth/sc/rna/mk

# Build genome
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--runThreadN 1 \
--runMode genomeGenerate \
--genomeDir ~/SCZ-BP/ \
--genomeFastaFiles ~/SCZ-BP/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz ~/SCZ-BP/gencode.vM23.primary_assembly.annotation.gtf.gz \

  
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/SCZ-BP/ \
--readFilesIn ${work_path2}/190620KelA_D19-6787_1_sequence.fastq.gz ${work_path2}/190620KelA_D19-6787_2_sequence.fastq.gz \
--soloType Droplet \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 16 \
--soloBarcodeReadLength 0 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures Gene
