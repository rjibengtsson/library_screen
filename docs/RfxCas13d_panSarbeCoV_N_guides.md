# Workflow for selection of RfxCas13 guides close to previously tested SarbeCoV N guides

Identify RfxCas13b guides predicted to be highly efficient by TIGER that overlap with a previously validated PspCas13b guide targeting pan-sarbecov nucleocapsid, which demonstrated >3 log fold reduction against VIC01 and Omicron variants.

First, identify where the 5 previously tested PspCas13b pan-Sarbecovirus guides map to the MERS N gene. These guides should be highly conserved among SarbeCoVs but will be divergent to MERS. Thus, blastn will have to allow mismatches. 

```bash
python3 data/raw/N/scripts/blastn.py
```

Grab RfxCas13d guides belonging to N region.

```bash
awk '$1 ~ /N/' data/raw/NC_019843.3/RfxCas13d_guides.bed > data/raw/RfxCas13d_guides.bed
```

PspCas13b bed file was manually generated from blastn output.

Generate RfxCas13d guides that are proximate to PspCas13b guide targeting SarbeCoV N gene.

```bash
python3 src/bedtools_closest.py -a data/raw/N/PspCas13b_guides.bed -b data/raw/N/RfxCas13d_guides.bed -o results -f Pan-SarbeCoV_guides.html
```

Once guides have been manually selected, extract guide sequence for RfxCas13d

```bash
python3 scripts/get_guide_seq.py -csv data/raw/N/closest_guides.csv -in data/raw/N/rfxcas13d_idx.txt -o data/raw/N/RfxCas13d_MERS_guides.csv
```

## SARS-CoV-2 nucleocapsid guides

The same pipeline was applied to select RfxCas13d guides over lapping with a PspCas13b guide targeting SARS-CoV-2 nucleocapsid, which demonstrated >3 log fold reduction against VIC01 and Omicron variants.
