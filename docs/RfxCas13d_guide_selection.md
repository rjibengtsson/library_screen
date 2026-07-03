# RfxCas13d guide selection workflow

## MERS

A MERS_PspCas13b_guides.csv file was provided by Khoa which contain guide seqeunces tiled across

Download MERS genome sequences and extract subgenomic regions. Modify lines 40 and 46 to specify your target subgenomic regions. 

The script will download GenBank and FASTA files for user-provided accessions from NCBI, extracts subgenomic region coordinates from the GenBank file, and retrieves the corresponding nucleotide sequences from the FASTA file.
```bash
python scripts/extract_subgenomic_seq.py -a NC_019843.3 -o data/raw/
``` 

Run BLASTN to map PspCas13b and RfxCas13d guides against a reference MERS genome and identify their positions.

This pipeline takes sequences from previously tested PspCas13b and RfxCas13d guides and maps them to a reference genome. It outputs two BED files to identify guides that intersect. 
```bash
python3 scripts/get_guide_positions.py -b data/MERS_PspCas13b_guides.csv -d data/mers_RfxCas13d_TIGER.csv -a NC_019843.3 -r data/raw/ -o data/raw/
```

Filter RxfCas13d bed file based on guide score and intersect with PspCas13b bed file to get intersection.

This script outputs a final BED file containing, for each previously tested PspCas13b guide, the closest RfxCas13d guide with the highest prediction score. 
```bash
python3 src/bedtools_closest.py -b data/processed/PspCas13b_guides.bed -a data/processed/RfxCas13d_guides.bed -o data/processed
qstart
```

## SARS-CoV-2

The same pipeline was used to identify guides predicted to be most potent for RfxCas13d targeting SARS-CoV-2, and closest in position to previously tested PspCas13b guides.

Download SARS-CoV-2 genome sequences. This time only  modify lines 40 and 46 to specify your target subgenomic regions. 

```bash
python scripts/extract_subgenomic_seq.py -a MT007544.1 -o data/raw/
```

Run BLASTN to map PspCas13b and RfxCas13d guides against a reference SARS-CoV-2 genome and identify their positions.

```bash
python3 scripts/get_guide_positions.py  -b data/SARS-CoV-2_PspCas13b_guides.xlsx -d data/sars_cov2_RfxCas13d_TIGER.xlsx  -r data/raw/MT007544.1/MT007544.1.fasta -o data/raw/
```

Filter RxfCas13d bed file based on guide score and intersect with PspCas13b bed file to get intersection.

```bash
python3 src/bedtools_closest.py -a data/raw/MT007544.1/PspCas13b_guides.bed -b data/raw/MT007544.1/RfxCas13d_guides.bed -o data/results -f SARS-CoV-2_guides.html
```
