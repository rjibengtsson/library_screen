"""
This script extracts subgenomic sequences from a given genome based on the provided GenBank accession number. 
It fetches the genome in both GenBank and FASTA formats, identifies feature locations, and extracts the corresponding sequences.
Author: Rebecca Bengtsson
Date: 2024-06-20
"""

from pyexpat import features
from Bio import SeqIO
from pathlib import Path
import sys
import argparse


# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(target_dir))
import src.entrez as entrez

def parse_args():
    parser = argparse.ArgumentParser(description="Extract subgenomic sequences from given genome")
    parser.add_argument("-a", "--accession", help="GenBank accession number", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Directory to save extracted sequences", type=Path, required=True)
    
    return parser.parse_args()




def fetch_files(accession, gbk_file_path, fasta_file_path):
    record = entrez.fetch_gbk(accession)
    SeqIO.write(record, gbk_file_path, "genbank")
    record = entrez.fetch_fasta(accession)
    SeqIO.write(record, fasta_file_path, "fasta")


def main():
    args = parse_args()
    accession = args.accession

    output_dir = args.output_dir / accession
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    gbk_file_path = output_dir / f"{accession}.gbk"
    fasta_file_path = output_dir / f"{accession}.fasta"

    # #### Step 1 - Download MERS genomes from NCBI
    # fetch_files(accession, gbk_file_path, fasta_file_path)

    feature_locations = entrez.get_feature_locations(gbk_file_path)
    
    #### Step 2 - Extract feature coordinates and save to file
    feature_list = [feature for feature in feature_locations.keys()]
    for feature in feature_list:
        start = feature_locations[feature]['start']
        end = feature_locations[feature]['end']
        entrez.fetch_gene_sequence(feature, start, end, fasta_file_path, output_dir)



if __name__ == "__main__":
    main()