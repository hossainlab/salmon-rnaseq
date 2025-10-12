#!/usr/bin/env bash

# Step 1. Download the latest GENCODE files
# Visit: https://www.gencodegenes.org/human/
# Fasta files > Transcript sequences > Right Click on Fasta > Copy link address > Paste here 
TRANS_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.pc_transcripts.fa.gz"


# Visit: https://www.gencodegenes.org/human/
# Fasta files > Genome sequence, primary assembly (GRCh38) > Right Click on Fasta > Copy link address > Paste here 
GENOME_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/GRCh38.primary_assembly.genome.fa.gz"

wget -c "${TRANS_URL}" -O ../inputs/gencode.v49.pc_transcripts.fa.gz
wget -c "${GENOME_URL}" -O ../inputs/GRCh38.primary_assembly.genome.fa.gz

# Step 2. Create the decoy list
grep '^>' <(gunzip -c ../inputs/GRCh38.primary_assembly.genome.fa.gz) | cut -d ' ' -f 1 > ../inputs/decoys.txt
sed -i -e 's/>//g' ../inputs/decoys.txt

# Step 3. Combine transcriptome and genome FASTAs
cat ../inputs/gencode.v49.pc_transcripts.fa.gz ../inputs/GRCh38.primary_assembly.genome.fa.gz > ../inputs/transcripts_and_decoys.fa.gz

# Step 4. Build the Salmon index
salmon index \
  -t ../inputs/transcripts_and_decoys.fa.gz \
  -d ../inputs/decoys.txt \
  -p 30 \
  -i ../inputs/human_salmon_index \
  --gencode

echo "All files and outputs saved in the existing ./input directory"
