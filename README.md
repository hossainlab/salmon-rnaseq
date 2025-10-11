# End-to-end bulk RNA-seq Quantification Pipeline using Salmon 

## Overview

This repository contains a **publication-ready RNA-seq quantification pipeline** built around **[Salmon](https://combine-lab.github.io/salmon/)** for lightweight and accurate transcript quantification.  
It automates the full workflow from **SRA download ? FASTQ generation ? Transcript quantification (Salmon)**, followed by **R-based downstream analysis** (DESeq2, functional enrichment, and visualization). The pipeline is designed for **reproducibility**, **modularity**, and **high-performance computing environments**, with detailed timing logs for benchmarking.

![](img/de_workflow_salmon.png)

## Features

- Fully automated and parallelized
- Ready for downstream DESeq2 analysis in R
- Conda-based reproducible environment via `environment.yml`
- Organized directory structure for clarity and reproducibility 



## Quick Start

```bash
conda env create -f environment.yml
conda activate salmon-rnaseq
```

