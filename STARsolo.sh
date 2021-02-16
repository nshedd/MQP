
work_path1=/data/zusers/sheddn
work_path2=/data/zusers/pratth/sc/rna/mk


# Then build index by STAR
sh test_build_index.sh 8 \
  ${work_path1}/SCZ-BP/mm10-2020-A_build/Mus_musculus.GRCm38.dna.primary_assembly.fa.modified \
# ${work_path1}/mm10-2020-A_build/gencode.vM23.primary_assembly.annotation.gtf.filtered \
  ${work_path1}/STAR_index_2.7.7a/
  
# Run STARsolo 2.7.7a
/data/zusers/pratth/sc/STAR/bin/Linux_x86_64/STAR \
--genomeDir ${work_path1}/STAR_index_2.7.7a/ \
--readFilesIn ${work_path2}/190620KelA_D19-6787_1_sequence.fastq.gz ${work_path2}/190620KelA_D19-6787_2_sequence.fastq.gz \
--soloType Droplet \
--soloCBwhitelist ${work_path1}/737K-august-2016.txt \
--runThreadN 16 \
--soloBarcodeReadLength 0 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes NH HI nM AS CR UR CB UB GX GN \
--soloFeatures Gene
