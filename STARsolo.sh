
work_path1=/data/zusers/sheddn
work_path2=/data/zusers/pratth/sc/rna/mk

genome="GRCh38"
version="2020-A"

build="GRCh38-2020-A_build"
mkdir -p "$build"

fasta_in="~/SCZ-BP/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf_in="~/SCZ-BP/gencode.v32.primary_assembly.annotation.gtf.gz"

fasta_modified="$build/$(basename "$fasta_in").modified"

cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"

# Build genome
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
--runThreadN 1 \
--runMode genomeGenerate \
--genomeDir ~/SCZ-BP/ \
--genomeFastaFiles ~/SCZ-BP/Homo_sapiens.GRCm38.dna.primary_assembly.fa.gz.modified
--sjdbGTFfile ~/SCZ-BP/gencode.v32.primary_assembly.annotation.gtf.gz.modified \

  
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
