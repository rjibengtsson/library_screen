# Cas13 library screen

This repository contains scripts and data for designing and selecting CRISPR-Cas13 gRNAs for a library screen targeting MERS-CoV and SARS-CoV-2 genome.

## Overview

The workflow supports guide RNA selection for two Cas13 variants — **PspCas13b** and **RfxCas13d** — against the MERS-CoV(NC_019843.3) and SARS-CoV-2 reference genome. 

The pipeline:

1. Downloads the MERS-CoV genome and annotation from NCBI and extracts subgenomic region sequences (ORFs, structural genes, UTRs).
2. Maps existing PspCas13b and RfxCas13d guide sequences to the reference genome using BLASTN and generates BED files of guide positions.
3. Identifies RfxCas13d guides with the highest prediction scores that are closest to previously validated PspCas13b guides, to prioritise candidates for a paired library screen.

## Repository structure

```
scripts/          # Command-line scripts for data processing
  extract_subgenomic_seq.py   # Download genome files and extract subgenomic sequences
  get_guide_positions.py      # Map guide RNAs to reference genome via BLASTN

src/              # Shared utility modules
  bedtools_closest.py         # Find closest high-scoring RfxCas13d guide per PspCas13b guide
  entrez.py                   # NCBI Entrez fetch helpers
  utils.py                    # General utilities

data/
  raw/            # Raw input files (genome FASTA/GenBank, BLAST DB, guide sequences)
  processed/      # Pipeline outputs (BED files, closest-guide table)

docs/             # Workflow documentation
notebooks/        # Exploratory analysis notebooks
```

## Usage

See [docs/RfxCas13d_guide_selection.md](docs/RfxCas13d_guide_selection.md) for the full step-by-step workflow.

## Requirements

Dependencies are listed in `environment.yaml`. To recreate the conda environment:

```bash
conda env create -f environment.yaml
```