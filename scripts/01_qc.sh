#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Input: FASTQ files in ../fastq/
# Output: QC reports in ../outputs/qc/

mkdir -p ../outputs/qc

echo "Running FastQC..."
fastqc ../fastq/*.fastq.gz -o ../outputs/qc/

echo "Running MultiQC..."
multiqc ../outputs/qc/ -o ../outputs/qc/

echo "Quality check complete! Reports saved in ../outputs/qc/"
