# Workflow for selection of RfxCas13 guides close to previously tested SarbeCoV N guides

Identify RfxCas13b guides predicted to be highly efficient (TIGER prediction) and also overlapping with 5 previously tested PspCas13b guides, designed to target SarbeCoV nucleocapsid protein.

First, identify where the 5 previously tested PspCas13b pan-Sarbecovirus guides map to the MERS N gene. These guides should be highly conserved among SarbeCoVs but will be divergent to MERS. Thus, blastn will have to allow mismatches. 

```bash
python3 data/raw/N/scripts/blastn.py
```

Filter for best PspCas13b guide hits (selecting highest bitscore)

```bash
python3 src/bedtools_closest.py -a data/raw/MT007544.1/PspCas13b_guides.bed -b data/raw/MT007544.1/RfxCas13d_guides.bed -o results -f SARS-CoV-2_guides.html
```