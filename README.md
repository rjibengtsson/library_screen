# Cas13 library screen

This repository contains scripts and data for designing and selecting CRISPR-Cas13 gRNAs for a library screen targeting MERS-CoV and SARS-CoV-2.

## Overview

The workflow supports guide RNA selection for two Cas13 variants — **PspCas13b** and **RfxCas13d** — against the MERS-CoV (NC_019843.3), SARS-CoV-2 (MT007544.1), and pan-sarbecovirus nucleocapsid (N gene) reference sequences.

The pipeline:

1. Downloads genome sequences and annotations from NCBI and extracts subgenomic region sequences (ORFs, structural genes, UTRs).
2. Maps PspCas13b and RfxCas13d guide sequences to the reference genome using BLASTN and generates BED files of guide positions.
3. Filters RfxCas13d guides by TIGER prediction score and identifies, for each PspCas13b guide, the closest high-scoring RfxCas13d guide to prioritise candidates for a paired library screen.
4. Extracts final guide sequences and metadata for selected candidates.

## Repository structure

```
scripts/
  extract_subgenomic_seq.py   # Download genome files (FASTA/GenBank) and extract subgenomic sequences
  get_guide_positions.py      # Map PspCas13b and RfxCas13d guides to a reference genome via BLASTN
  get_guide_seq.py            # Extract guide sequences and metadata for selected candidates

src/                          # Shared utility modules
  bedtools_closest.py         # Filter RfxCas13d guides by score and find closest per PspCas13b guide
  entrez.py                   # NCBI Entrez fetch helpers
  utils.py                    # General utilities
  visualisation.py            # Guide position visualisation (Plotly)

data/
  raw/
    NC_019843.3/              # MERS-CoV reference genome, BLAST DB, and guide BED files
    MT007544.1/               # SARS-CoV-2 reference genome, BLAST DB, and guide BED files
    N/                        # Pan-sarbecovirus N gene sequences and guide mapping outputs
  processed/                  # Final pipeline outputs (BED files, closest-guide tables)

results/                      # Interactive HTML guide visualisations
docs/                         # Step-by-step workflow documentation
notebooks/                    # Exploratory analysis notebooks
```

## Usage

Detailed step-by-step workflows are in the `docs/` folder:

- [docs/RfxCas13d_guide_selection.md](docs/RfxCas13d_guide_selection.md) — MERS-CoV and SARS-CoV-2 guide selection
- [docs/RfxCas13d_panSarbeCoV_N_guides.md](docs/RfxCas13d_panSarbeCoV_N_guides.md) — Pan-sarbecovirus N gene guide selection

### Quick reference

```bash
# 1. Download genome and extract subgenomic sequences
python scripts/extract_subgenomic_seq.py -a NC_019843.3 -o data/raw/

# 2. Map guides to reference genome
python scripts/get_guide_positions.py \
    -b data/MERS_PspCas13b_guides.csv \
    -d data/MERS_RfxCas13d_TIGER.xlsx \
    -r data/raw/NC_019843.3/NC_019843.3.fasta \
    -o data/raw/

# 3. Find closest high-scoring RfxCas13d guide per PspCas13b guide
python src/bedtools_closest.py \
    -a data/raw/NC_019843.3/PspCas13b_guides.bed \
    -b data/raw/NC_019843.3/RfxCas13d_guides.bed \
    -o results -f MERS_guides_v2.html

# 4. Extract sequences for selected guides
python scripts/get_guide_seq.py \
    -csv data/raw/NC_019843.3/closest_guides.csv \
    -in data/processed/rfxcas13d_idx.txt \
    -o results/RfxCas13d_MERS_guides.csv
```

## Requirements

Dependencies are listed in `environment.yaml`. To recreate the conda environment:

```bash
conda env create -f environment.yaml
conda activate library_screen
```