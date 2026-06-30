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
    gbk_file_path = args.output_dir / f"{accession}.gbk"
    fasta_file_path = args.output_dir / f"{accession}.fasta"

    #### Step 1 - Download MERS genomes from NCBI
    fetch_files(accession, gbk_file_path, fasta_file_path)

    #### Step 2 - Extract feature coordinates and save to file
    # Get feature types
    feature_type_list = ["5'UTR", "3'UTR"]
    for feature in feature_type_list:
        results = entrez.fetch_gene_coordinates(gbk_file_path, feature_type=feature)
        entrez.save_record(results, "feature_type", fasta_file_path, args.output_dir)

    # Get genes
    gene_list = ["orf1ab", "S", "N", "E", "M"]
    for gene in gene_list:
        results = entrez.fetch_gene_coordinates(gbk_file_path, gene_name=gene)
        entrez.save_record(results, "gene", fasta_file_path, args.output_dir)


if __name__ == "__main__":
    main()