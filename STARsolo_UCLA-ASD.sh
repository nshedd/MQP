work_path1=/data/rusers/sheddn/UCLA-ASD/data
work_path2=/data/projects/psychencode/TR-SingleCell/UCLA-ASD/PEC_syn18898607/scRNAseq

# echo "10BW_S13"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/10BW_S13_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/10BW_S13_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/10BW_S13_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/10BW_S13_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/737K-august-2016.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/10BW_S13/
#   
#   echo "11BW_S14"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/11BW_S14_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/11BW_S14_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/11BW_S14_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/11BW_S14_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/737K-august-2016.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/11BW_S14/
#   
#   echo "12BW_S15"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/12BW_S15_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/12BW_S15_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/12BW_S15_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/12BW_S15_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/737K-august-2016.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/12BW_S15/
#   
#   echo "13BW_S16"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/13BW_S16_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/13BW_S16_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/13BW_S16_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/13BW_S16_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/737K-august-2016.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/13BW_S16/

# echo "14BW_S17"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/14BW_S17_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/14BW_S17_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/14BW_S17_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/14BW_S17_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/737K-august-2016.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/14BW_S17/
#   
#   echo "15BW_S1"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/15BW_S1_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/15BW_S1_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/15BW_S1_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/15BW_S1_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/3M-february-2018.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --soloUMIlen 12 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/15BW_S1/
#   
#   echo "16BW_S2"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/16BW_S2_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/16BW_S2_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/16BW_S2_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/16BW_S2_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/3M-february-2018.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --soloUMIlen 12 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/16BW_S2/
#   
#   echo "17BW_S3"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/17BW_S3_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/17BW_S3_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/17BW_S3_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/17BW_S3_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/3M-february-2018.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --soloUMIlen 12 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/17BW_S3/
#   
#   echo "18BW_S4"
# # Run STARsolo 2.7.7a
# ~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
# --genomeDir ~/UCLA-ASD/ \
# --readFilesIn <(gunzip -c ${work_path2}/18BW_S4_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/18BW_S4_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/18BW_S4_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/18BW_S4_L002_R1_001.fastq.gz) \
# --soloType CB_UMI_Simple \
# --soloCBwhitelist ~/3M-february-2018.txt \
# --runThreadN 12 \
# --soloBarcodeReadLength 1 \
# --soloUMIlen 12 \
# --outSAMtype BAM SortedByCoordinate \
# --outSAMunmapped Within \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
# --soloFeatures GeneFull \
# --outFileNamePrefix ${work_path1}/18BW_S4/
  
  echo "19BW_S5"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/19BW_S5_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/19BW_S5_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/19BW_S5_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/19BW_S5_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/19BW_S5/
  
  echo "1BW_S1"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/1BW_S1_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/1BW_S1_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/1BW_S1_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/1BW_S1_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/1BW_S1/
  
  
  echo "20BW_S6"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/20BW_S6_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/20BW_S6_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/20BW_S6_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/20BW_S6_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/20BW_S6/
  
  
  echo "21BW_S7"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/21BW_S7_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/21BW_S7_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/21BW_S7_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/21BW_S7_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/21BW_S7/
  
  
  echo "22BW_S8"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/22BW_S8_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/22BW_S8_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/22BW_S8_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/22BW_S8_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/22BW_S8/
  
  
  echo "23BW_S9"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/23BW_S9_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/23BW_S9_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/23BW_S9_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/23BW_S9_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/23BW_S9/
  
  echo "24BW_S10"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/24BW_S10_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/24BW_S10_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/24BW_S10_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/24BW_S10_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/24BW_S10/
  
  echo "25BW_S1"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/25BW_S1_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/25BW_S1_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/25BW_S1_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/25BW_S1_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/25BW_S1/
  
  echo "26BW_S2"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/26BW_S2_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/26BW_S2_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/26BW_S2_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/26BW_S2_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/26BW_S2/
  
  echo "27BW_S3"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/27BW_S3_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/27BW_S3_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/27BW_S3_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/27BW_S3_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/27BW_S3/
  
  echo "28BW_S4"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/28BW_S4_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/28BW_S4_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/28BW_S4_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/28BW_S4_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/28BW_S4/
  
  echo "29BW_S5"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/29BW_S5_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/29BW_S5_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/29BW_S5_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/29BW_S5_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/29BW_S5/
  
  echo "2BW_S2"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/2BW_S2_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/2BW_S2_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/2BW_S2_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/2BW_S2_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/2BW_S2/
  
  echo "30BW_S6"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/30BW_S6_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/30BW_S6_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/30BW_S6_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/30BW_S6_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/30BW_S6/
  
  echo "31BW_S1"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/31BW_S1_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/31BW_S1_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/31BW_S1_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/31BW_S1/
  
  echo "32BW_S2"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/32BW_S2_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/32BW_S2_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/32BW_S2_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/32BW_S2/
  
  echo "33BW_S3"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/33BW_S3_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/33BW_S3_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/33BW_S3_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/33BW_S3/
  
  echo "34BW_S4"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/34BW_S4_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/34BW_S4_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/34BW_S4_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/34BW_S4/
  
  echo "35BW_S5"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/35BW_S5_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/35BW_S5_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/35BW_S5_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/35BW_S5/
  
  echo "37BW_S7"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/37BW_S7_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/37BW_S7_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/37BW_S7_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/37BW_S7/
  
  echo "38BW_S8" 
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/38BW_S8_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/38BW_S8_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/38BW_S8_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/38BW_S8/
  
  echo "39BW_S9"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/39BW_S9_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/39BW_S9_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/39BW_S9_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/39BW_S9/
  
  echo "3BW_S3"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/3BW_S3_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/3BW_S3_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/3BW_S3_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/3BW_S3_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/3BW_S3/
  
  echo "40BW_S10"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/40BW_S10_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/40BW_S10_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/40BW_S10_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/40BW_S10/
  
  echo "41BW_S11"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/41BW_S11_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/41BW_S11_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/41BW_S11_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/41BW_S11/
  
  echo "42BW_S12"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/42BW_S12_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/42BW_S12_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/42BW_S12_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/42BW_S12/
  
  echo "43BW_S13"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/43BW_S13_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/43BW_S13_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/43BW_S13_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/43BW_S13/
  
  echo "44BW_S14"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/44BW_S14_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/44BW_S14_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/44BW_S14_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/44BW_S14/
  
  echo "45BW_S15"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/45BW_S15_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/45BW_S15_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/45BW_S15_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/45BW_S15/
  
  echo "46BW_S1"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/46BW_S1_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/46BW_S1_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/46BW_S1_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/46BW_S1/
  
  echo "47BW_S2"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/47BW_S2_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/47BW_S2_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/47BW_S2_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/47BW_S2/
  
  echo "48BW_S3"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/48BW_S3_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/48BW_S3_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/48BW_S3_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/48BW_S3/
  
  echo "49BW_S4"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/49BW_S4_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/49BW_S4_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/49BW_S4_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/49BW_S4/
  
  echo "4BW_S4"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/4BW_S4_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/4BW_S4_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/4BW_S4_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/4BW_S4_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/4BW_S4/
  
  echo "50BW_S5"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/50BW_S5_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/50BW_S5_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/50BW_S5_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/50BW_S5/
  
  echo "51BW_S6"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/51BW_S6_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/51BW_S6_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/51BW_S6_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/51BW_S6/
  
  echo "52BW_S7"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/52BW_S7_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/52BW_S7_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/52BW_S7_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/52BW_S7/
  
  echo "53BW_S8"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/53BW_S8_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/53BW_S8_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/53BW_S8_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/53BW_S8/
  
  echo "54BW_S11"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/54BW_S11_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/54BW_S11_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/54BW_S11_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/54BW_S11/
  
  echo "55BW_S12"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/55BW_S12_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/55BW_S12_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/55BW_S12_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/55BW_S12/
  
  echo "56BW_S13"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/56BW_S13_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/56BW_S13_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/56BW_S13_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/56BW_S13/
  
  echo "57BW_S14"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/57BW_S14_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/57BW_S14_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/57BW_S14_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/57BW_S14/
  
  echo "58BW_S28"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/58BW_S28_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/58BW_S28_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/58BW_S28_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/58BW_S28/
  
  echo "59BW_S29"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/59BW_S29_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L002_R2_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L003_R2_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L004_R2_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L005_R2_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L006_R2_001.fastq.gz) <(gunzip -c ${work_path2}/59BW_S29_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L002_R1_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L003_R1_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L004_R1_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L005_R1_001.fastq.gz),<(gunzip -c ${work_path2}/59BW_S29_L006_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/3M-february-2018.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--soloUMIlen 12 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/59BW_S29/
  
  echo "5BW_S5"
# Run STARsolo 2.7.7a
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--genomeDir ~/UCLA-ASD/ \
--readFilesIn <(gunzip -c ${work_path2}/5BW_S5_L001_R2_001.fastq.gz),<(gunzip -c ${work_path2}/5BW_S5_L002_R2_001.fastq.gz) <(gunzip -c ${work_path2}/5BW_S5_L001_R1_001.fastq.gz),<(gunzip -c ${work_path2}/5BW_S5_L002_R1_001.fastq.gz) \
--soloType CB_UMI_Simple \
--soloCBwhitelist ~/737K-august-2016.txt \
--runThreadN 12 \
--soloBarcodeReadLength 1 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures GeneFull \
--outFileNamePrefix ${work_path1}/5BW_S5/
  
  
  
