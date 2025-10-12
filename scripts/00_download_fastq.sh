#!/usr/bin/env bash

# Download FASTQ files using fastq-dl
# Reads accessions from ../inputs/SRR_Acc_List.txt
# Saves all FASTQ files in the root-level fastq/ directory

# Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

mkdir -p ../fastq

while read -r ACC; do
    [[ -z "$ACC" ]] && continue
    echo "Downloading $ACC ..."
    fastq-dl --accession "$ACC" --outdir ../fastq
done < ../inputs/SRR_Acc_List.txt

echo "All FASTQ files are saved in ../fastq/"
