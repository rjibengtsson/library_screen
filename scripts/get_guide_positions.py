"""
This script identifies guide RNA positions for Cas13b and Cas13d systems by BLASTing guide 
sequences against a reference FASTA, then outputs BED files with the matched positions.
Author: Rebecca Bengtsson
Date: 2024-06-20
"""


import os
import sys
from pathlib import Path
import argparse

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(target_dir))
import src.utils as utils
import src.clean_resource as clean_resource



def parse_args():
    parser = argparse.ArgumentParser(
        description="Perform BLASTN to search for guide RNA positions"
    )
    
    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-b', '--cas13b', 
        type=str, required=True,
        help='Cas13b guide sequence file (accept [*.csv, *.tsv, *.xlsx] format)'
    )
    required.add_argument(
        '-d', '--cas13d',
        type=str, required=True,
        help='Cas13d guide sequence file (accept [*.csv, *.tsv, *.xlsx] format)'
    )
    required.add_argument(
        '-o', '--output_dir',
        type=Path, required=True,
        help='Output directory for BLASTN results'
    )
    
    # Input source (mutually exclusive)
    input_source = parser.add_argument_group('input source', 'Choose one:')
    input_group = input_source.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-r', '--reference',
        type=str,
        help='Reference genome FASTA file'
    )
    input_group.add_argument(
        '-db', '--database',
        type=str,
        help='BLASTN database FASTA file'
    )
    
    # Optional arguments
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument(
        '-p', '--pident',
        type=int, default=100,
        help='Percentage identity threshold for BLASTN (default: 100)'
    )
    optional.add_argument(
        '-c', '--cov_perc',
        type=int, default=100,
        help='Percentage coverage threshold for BLASTN (default: 100)'
    )
    
    return parser.parse_args()




def _extract_accession(filename: str) -> str:
    """Extract accession ID by removing common sequence file extensions."""
    return filename.removesuffix('.fasta').removesuffix('.fa').removesuffix('.fna')



def main():
    args = parse_args()
    cas13b_file = args.cas13b
    cas13d_file = args.cas13d
    
    # Step 1 - Determine BLAST database and output directory
    if args.database is not None:
        db_path = Path(str(Path(args.database)).removesuffix('.fasta').removesuffix('.fa').removesuffix('.fna'))
        output_dir = db_path.parent
        accession = _extract_accession(db_path.name)
        # db_path = output_dir / f"{accession}_subgenomic.fasta"
        print(f"Using existing BLAST database: {db_path}")
    
    elif args.reference is not None:
        reference_path = Path(args.reference)
        output_dir = reference_path.parent
        db_path = Path(utils.make_blast_db(args.reference))
        accession = _extract_accession(reference_path.name)
    
    else:
        raise ValueError("Either --reference or --database must be specified.")

    # Step 2 - Standardise and combine guide data
    cas13b_df = clean_resource.clean_cas13b(cas13b_file)
    cas13d_df = clean_resource.clean_cas13d(cas13d_file)
    combined_df = clean_resource.combine_guides(cas13b_df, cas13d_df, output_dir)

    # Step 3 - Generate guide FASTA
    query_fasta = utils.generate_guide_fasta(combined_df, output_dir)


    # set default values for pident and cov_perc if not provided
    if args.pident is not None:
        pident = args.pident
    else:
        pident = 100

    if args.cov_perc is not None:
        cov_perc = args.cov_perc
    else:
        cov_perc = 100


    # Step 4 - Run BLASTN
    blastn_output_file = output_dir / "blastn_output.txt"
    utils.blastn_query(query_fasta, db_path, blastn_output_file, pident=pident, cov_perc=cov_perc)

    # Step 5 - Process BLASTN results and generate BED files
    utils.generate_bed_file(blastn_output_file, combined_df, output_dir)



if __name__ == "__main__":
    main()