"""
This module contains utility functions for processing genomic data, 
including creating BLAST databases, performing BLAST searches, 
and generating output files in various formats.
"""

import subprocess
import os, sys
import glob
from pathlib import Path
from tracemalloc import start
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
import entrez as entrez



def make_blast_db(reference_dir):

    base = os.path.basename(reference_dir)
    accession = base.removesuffix('.fasta').removesuffix('.fa').removesuffix('.fna')

    ref_dir = os.path.dirname(reference_dir)

    # Find all *.fasta files, excluding those starting with the accession number
    pattern = os.path.join(ref_dir, "*.fasta")
    fasta_files = sorted(
        f for f in glob.glob(pattern)
        if not os.path.basename(f).startswith(accession)
    )

    concat_path = os.path.join(ref_dir, f"{accession}_subgenomic.fasta")

    # Concatenate all matching FASTA files
    with open(concat_path, "w") as outfile:
        for fasta_file in fasta_files:
            with open(fasta_file, "r") as infile:
                outfile.write(infile.read())

    print(f"Concatenated to: {concat_path}")

    # Build BLAST nucleotide database
    db_path = os.path.splitext(concat_path)[0]
    subprocess.run(
        ["makeblastdb", "-in", concat_path, "-dbtype", "nucl", "-out", db_path],
        check=True,
    )
    return db_path


def blastn_query(query_fasta, db_path, output_file, pident=100, cov_perc=100):
    blastn_cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", f"{query_fasta}",
        "-db", f"{db_path}",
        "-max_target_seqs", "100",
        "-perc_identity", f"{pident}",
        "-qcov_hsp_perc", f"{cov_perc}",
        "-outfmt", "6",
        "-out", f"{output_file}"]
    print(f"Running BLASTN with command: {' '.join(blastn_cmd)}")
    subprocess.run(blastn_cmd, check=True)



def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]



def file_to_df(file_path):
    file_path = Path(file_path)
    # Detect the file extension and read accordingly
    if file_path.suffix == '.csv':
        df = pd.read_csv(file_path, header=0)
    elif file_path.suffix in ['.tsv', '.txt']:
        df = pd.read_csv(file_path, sep='\t', header=0)
    elif file_path.suffix in ['.xls', '.xlsx']:
        sheets = pd.read_excel(file_path, sheet_name=None, header=0) # Read all the sheets
        df = pd.concat(sheets.values(), ignore_index=True)
        # # print(pd.read_excel(file_path,header=0))
        # print(df)
    else:
        raise ValueError("Unsupported file format.")
    return df   


def generate_guide_fasta(df, output_dir):
    records = [
        SeqRecord(Seq(seq), id=str(seq_id), description="")
        for seq_id, seq in zip(df['Cas'] + '_' + df['Index'].astype(str), df["Target Sequence"])
    ]

    SeqIO.write(records, output_dir / "query_guides.fasta", "fasta")
    return output_dir / "query_guides.fasta"


def filter_blastn_results(blastn_output):
    # for each sseqid, keep the row with the highest bitscore
    df = pd.read_csv(blastn_output, sep="\t", header=None,
                            names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                   "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    df = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
    # df.to_csv(blastn_output, sep="\t", header=False, index=False)
    return df




def generate_bed_file(blastn_output, combined_df, output_dir):
    blastn_df = filter_blastn_results(blastn_output)
    # print(blastn_df[blastn_df['qseqid'].str.startswith('PspCas13b')])
    
    # Map Guide sequence and target sequence to the blastn_df
    blastn_df[['Cas', 'Index']] = blastn_df['qseqid'].str.split('_', n=1, expand=True)
    blastn_df['Index'] = blastn_df['Index'].astype(int)
    blastn_df['gene'] = blastn_df['sseqid'].str.split('_').str[-1]
    blastn_df['sstart'] -= 1  # Convert to 0-based

    # Merge Guide Sequence and Target Sequence on Cas + Index
    seq_cols = combined_df[['Cas', 'Index', 'Guide Sequence', 'Target Sequence']]
    blastn_df = blastn_df.merge(seq_cols, on=['Cas', 'Index'], how='left')

    # Map RfxCas13d guide scores by index
    scores = combined_df[combined_df['Cas'] == 'RfxCas13d'].set_index('Index')['Guide Score']
    rfx_mask = blastn_df['Cas'] == 'RfxCas13d'
    blastn_df.loc[rfx_mask, 'Guide Score'] = blastn_df.loc[rfx_mask, 'Index'].map(scores)

    # Write BED files per Cas type
    for cas, group in blastn_df.groupby('Cas'):
        cols = ['sseqid', 'sstart', 'send', 'Index', 'Guide Sequence', 'Target Sequence']
        if cas == 'RfxCas13d':
            cols.append('Guide Score')
        group[cols].to_csv(output_dir / f"{cas}_guides.bed", sep='\t', header=False, index=False)
