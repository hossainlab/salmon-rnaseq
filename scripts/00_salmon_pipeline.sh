#!/bin/bash

# --------------------------------------------------------------
# Author: Md. Jubayer Hossain
# Organization: DeepBio Limited / CHIRAL Bangladesh
# Project: RNA-seq Quantification with Salmon
# Description: SRA ? FASTQ ? Salmon ? Quantification
# Date: 2025-10-06
# Version: 1.2 (repo-aware paths)

# --------------------------------------------------------------
# Data 
# GEO link: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778
# SRA link: https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=229998

samples_to_quant=(SRR1039508 SRR1039509 SRR1039512 SRR1039513 SRR1039516 SRR1039517 SRR1039520)

for sample in "${samples_to_quant[@]}"
do
  
  # Print the current sample
  echo $sample
  
  # Prefetch
  prefetch -O $sample/ $sample
  
  # Fasterq-Dump
  fasterq-dump -e 14 -p -O $sample/ $sample/$sample/$sample.sra
  
  # Salmon quant
  salmon quant -l A -1 $sample/$sample.sra_1.fastq -2 $sample/$sample.sra_2.fastq  --validateMappings -i ~/hg38/salmon_partial_sa_index/default/ -o Salmon.out/$sample -p 15
  
done