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
    parser = argparse.ArgumentParser(description="Perform blastn to search for guide RNA positions")
    parser.add_argument("-b", "--cas13b", help="CAS13b sequence file", type=str, required=True)
    parser.add_argument("-d", "--cas13d", help="CAS13d sequence file", type=str, required=True)
    parser.add_argument("-r", "--reference", help="Reference sequence fasta file", type=str, required=True)
    parser.add_argument("-o", "--output_dir", help="Directory to save blastn output", type=Path, required=True)
    
    try:
        args = parser.parse_args()
        if not args.cas13b or not args.cas13d or not args.reference or not args.output_dir:
            parser.print_help()
            sys.exit(1)
    except Exception as e:
        print(f"ERROR:Please check path exists: {e}")
        parser.print_help()
        sys.exit(1)

    return parser.parse_args()


def main():
    args = parse_args()
    cas13b_file = args.cas13b
    cas13d_file = args.cas13d
    reference_dir = args.reference
    output_dir = Path(os.path.dirname(reference_dir))

    #### Step 1 - generate blast database from reference sequences
    db_path = utils.make_blast_db(reference_dir)

    #### Step 2.1 - standardise Cas13b data
    cas13b_df = clean_resource.clean_cas13b(cas13b_file)

    #### Step 2.2 - standardise Cas13d data
    cas13d_df = clean_resource.clean_cas13d(cas13d_file)


    #### Step 2.3 - combine Cas13b and Cas13d data
    combined_df = clean_resource.combine_guides(cas13b_df, cas13d_df, output_dir)

    ### Step 3 - perform blastn for guide sequences and save output
    query_fasta = utils.generate_guide_fasta(combined_df, output_dir)


    #### Step 4 - perform blastn for guide sequences and save output
    accession = os.path.basename(reference_dir).removesuffix('.fasta').removesuffix('.fa').removesuffix('.fna')
    blastn_output_file = output_dir / f"blastn_output.txt"
    utils.blastn_query(query_fasta, db_path, blastn_output_file)

    #### Step 4 - make bed files for PspCas13b and RfxCas13d guides
    utils.generate_bed_file(blastn_output_file, combined_df, output_dir)




if __name__ == "__main__":
    main()