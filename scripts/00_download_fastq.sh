#!/usr/bin/env bash

# Download FASTQ files using fastq-dl
# Reads accessions from ../inputs/SRR_Acc_List.txt
# Saves all FASTQ files in the root-level fastq/ directory
# Example:
# parallel-fastq-dump --sra-id SRR2244401 --threads 4 --outdir out/ --split-files --gzip

# Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

mkdir -p ../fastq

while read -r ACC; do
    [[ -z "$ACC" ]] && continue
    echo "Downloading $ACC ..."
    parallel-fastq-dump \
        --sra-id "$ACC" \
        --threads 4 \
        --outdir ../fastq \
        --split-files \
        --gzip
done < ../inputs/SRR_Acc_List.txt

echo "All FASTQ files are saved in ../fastq/"