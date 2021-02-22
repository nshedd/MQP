work_path1=/home/sheddn

genome="GRCh38"
version="2020-A"

build="GRCh38-2020-A_build"
mkdir -p "$build"

fasta_in="/home/sheddn/UCLA-ASD/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
gtf_in="/home/sheddn/UCLA-ASD/gencode.v29.primary_assembly.annotation_UCSC_names.gtf"

fasta_modified="$build/$(basename "$fasta_in").modified"

ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
 cat "$fasta_in" \
     | sed -E 's/^>(\S+).*/>\1 \1/' \
     | sed -E 's/^>([0-9]+|[XY]) />chr\1 /' \
     | sed -E 's/^>MT />chrM /' \
     > "$fasta_modified"
     
# Remove version suffix from transcript, gene, and exon IDs in order to match
# previous Cell Ranger reference packages
#
# Input GTF:
#     ... gene_id "ENSG00000223972.5"; ...
# Output GTF:
#     ... gene_id "ENSG00000223972"; gene_version "5"; ...
gtf_modified="$build/$(basename "$gtf_in").modified"
# Pattern matches Ensembl gene, transcript, and exon IDs for human or mouse:
ID="(ENS(MUS)?[GTE][0-9]+)\.([0-9]+)"
cat "$gtf_in" \
    | sed -E 's/gene_id "'"$ID"'";/gene_id "\1"; gene_version "\3";/' \
    | sed -E 's/transcript_id "'"$ID"'";/transcript_id "\1"; transcript_version "\3";/' \
    | sed -E 's/exon_id "'"$ID"'";/exon_id "\1"; exon_version "\3";/' \
    > "$gtf_modified"

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

 # Create reference package
 /home/sheddn/yard/apps/cellranger-5.0.1/cellranger mkref --ref-version="$version" \
     --genome="$genome" --fasta="$fasta_modified"

# Build genome
~/STAR-2.7.7a/bin/Linux_x86_64/STAR \
 --runThreadN 1 \
 --runMode genomeGenerate \
 --genomeDir ~/UCLA-ASD \
 --genomeFastaFiles /home/sheddn/git/MQP/GRCh38-2020-A_build/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.modified \
 --sjdbGTFfile /home/sheddn/git/MQP/GRCh38-2020-A_build/gencode.v29.primary_assembly.annotation_UCSC_names.gtf.modified \
 --outFileNamePrefix /home/sheddn/UCLA-ASD \
