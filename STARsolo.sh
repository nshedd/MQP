
work_path1=/data/zusers/sheddn
work_path2=/data/zusers/pratth/sc/rna/mk

genome="GRCh38"
version="2020-A"

build="GRCh38-2020-A_build"
mkdir -p "$build"

fasta_in="/home/sheddn/SCZ-BP/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf_in="/home/sheddn/SCZ-BP/gencode.v32.primary_assembly.annotation.gtf.gz"

fasta_modified="$build/$(basename "$fasta_in").modified"

ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$fasta_in" \
    | sed -E 's/^>(\S+).*/>\1 \1/' \
    | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
    | sed -E 's/^>MT />chrM /' \
    > "$fasta_modified"
    
# Define string patterns for GTF tags
# NOTES:
# - Since GENCODE release 31/M22 (Ensembl 97), the "lincRNA" and "antisense"
#   biotypes are part of a more generic "lncRNA" biotype.
# - These filters are relevant only to GTF files from GENCODE. The GTFs from
#   Ensembl release 98 have the following differences:
#   - The names "gene_biotype" and "transcript_biotype" are used instead of
#     "gene_type" and "transcript_type".
#   - Readthrough transcripts are present but are not marked with the
#     "readthrough_transcript" tag.
#   - Only the X chromosome versions of genes in the pseudoautosomal regions
#     are present, so there is no "PAR" tag.
BIOTYPE_PATTERN=\
"(protein_coding|lncRNA|\
IG_C_gene|IG_D_gene|IG_J_gene|IG_LV_gene|IG_V_gene|\
IG_V_pseudogene|IG_J_pseudogene|IG_C_pseudogene|\
TR_C_gene|TR_D_gene|TR_J_gene|TR_V_gene|\
TR_V_pseudogene|TR_J_pseudogene)"
GENE_PATTERN="gene_type \"${BIOTYPE_PATTERN}\""
TX_PATTERN="transcript_type \"${BIOTYPE_PATTERN}\""
READTHROUGH_PATTERN="tag \"readthrough_transcript\""
PAR_PATTERN="tag \"PAR\""


# Construct the gene ID allowlist. We filter the list of all transcripts
# based on these criteria:
#   - allowable gene_type (biotype)
#   - allowable transcript_type (biotype)
#   - no "PAR" tag (only present for Y chromosome PAR)
#   - no "readthrough_transcript" tag
# We then collect the list of gene IDs that have at least one associated
# transcript passing the filters.
cat "$gtf_modified" \
    | awk '$3 == "transcript"' \
    | grep -E "$GENE_PATTERN" \
    | grep -E "$TX_PATTERN" \
    | grep -Ev "$READTHROUGH_PATTERN" \
    | grep -Ev "$PAR_PATTERN" \
    | sed -E 's/.*(gene_id "[^"]+").*/\1/' \
    | sort \
    | uniq \
    > "${build}/gene_allowlist"


# Filter the GTF file based on the gene allowlist
gtf_filtered="${build}/$(basename "$gtf_in").filtered"
# Copy header lines beginning with "#"
grep -E "^#" "$gtf_modified" > "$gtf_filtered"
# Filter to the gene allowlist
grep -Ff "${build}/gene_allowlist" "$gtf_modified" \
    >> "$gtf_filtered"

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
