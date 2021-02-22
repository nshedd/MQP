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
