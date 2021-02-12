
work_path=/data/zusers/sheddn/Geschwind

# download data (this part is the same as cellranger)
# https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#mm10_#{files.refdata_mm10.version}

# Then build index by STAR
sh test_build_index.sh 8 \
  ${work_path}/mm10-2020-A_build/Mus_musculus.GRCm38.dna.primary_assembly.fa.modified \
  ${work_path}/mm10-2020-A_build/gencode.vM23.primary_assembly.annotation.gtf.filtered \
  ${work_path}/STAR_index_2.7.7a/
  
# Run STARsolo 2.7.7a
STAR \
--genomeDir ${work_path}/STAR_index_2.7.7a/ \
--readFilesIn ${work_path}/data/ENCSR874BOF_S1_L001_R2_001.fastq ${work_path}/data/ENCSR874BOF_S1_L001_R1_001.fastq \
--soloType Droplet \
--soloCBwhitelist /home/liying/cellranger-5.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt \
--runThreadN 16 \
--soloBarcodeReadLength 0 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures Gene
