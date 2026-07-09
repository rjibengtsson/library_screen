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


def blastn_query(query_fasta, db_path, output_file, pident=70, cov_perc=40):
    blastn_cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", f"{query_fasta}",
        "-db", f"{db_path}",
        "-max_target_seqs", "100",
        "-perc_identity", f"{pident}",
        "-qcov_hsp_perc", f"{cov_perc}",
        "-outfmt", "6",
        "-word_size", "6",
        "-gapopen", "0",
        "-gapextend", "2",
        "-reward", "1",
        "-penalty", "-1",
        "-out", f"{output_file}"]
    print(f"Running BLASTN with command: {' '.join(blastn_cmd)}")
    subprocess.run(blastn_cmd, check=True)



def filter_blastn_results(blastn_output, final_output_file):
    # for each sseqid, keep the row with the highest bitscore
    df = pd.read_csv(blastn_output, sep="\t", header=None,
                            names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                   "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    # df = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
    df.to_csv(final_output_file, sep="\t", index=False)



def main():
    # Example usage
    home_dir = Path(__file__).resolve().parent.parent
    query_fasta = home_dir / "PspCas13b_query.fasta"
    db_path = home_dir / "N"
    output_file = home_dir / "blastn_result.txt"
    blastn_query(query_fasta, db_path, output_file)
    final_output_file = home_dir / "SarsCoV2_N_blastn_nogap.txt"
    filter_blastn_results(output_file, final_output_file)


if __name__ == "__main__":
    main()