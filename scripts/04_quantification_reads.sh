#!/bin/bash

# Set the path to the Salmon index
salmon_index="../inputs/human_salmon_index"

# Set the path to the "fastq" folder
fastq_dir="../fastq"

# Loop through R1 files in the fastq folder
for r1_file in "${fastq_dir}"/SRR*_1.fastq.gz; do
    # Derive the R2 file and sample name
    samp=$(basename "${r1_file%_1.fastq.gz}")
    r2_file="${fastq_dir}/${samp}_2.fastq.gz"

    echo "Processing sample ${samp}"
    salmon quant -i "$salmon_index" -l A \
        -1 "$r1_file" \
        -2 "$r2_file" \
        -p 28 --validateMappings -o "../outputs/salmon_out/${samp}"
done
