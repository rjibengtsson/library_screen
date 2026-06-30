import sys
from pathlib import Path
import argparse
import pandas as pd

# Add the directory containing your other modules
target_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(target_dir))
import src.utils as utils


def parse_args():
    parser = argparse.ArgumentParser(description="Perform blastn to search for guide RNA positions")
    parser.add_argument("-b", "--cas13b", help="CAS13b sequence file", type=str, required=True)
    parser.add_argument("-d", "--cas13d", help="CAS13d sequence file", type=str, required=True)
    parser.add_argument("-a", "--accession", help="Reference accession number", type=str, required=True)
    parser.add_argument("-r", "--reference", help="Reference sequence file directory", type=str, required=True)
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
    output_dir = args.output_dir

    #### Step 1 - generate blast database from reference sequences
    # db_path = utils.make_blast_db(args.accession, reference_dir)

    #### Step 2.1 - standardise Cas13b data
    cas13b_df = utils.read_csv(cas13b_file)
    cas13b_df['Index'] = cas13b_df.index
    cas13b_df = cas13b_df.rename(columns={"Spacer": "Guide Sequence", "Gene": "Region"})

    replacement = {
            "ORF1a": "ORF1ab",
            "ORF1a/b": "ORF1ab",
            "ORF1b": "ORF1ab"
        }

    cas13b_df['Region'] = cas13b_df['Region'].replace(replacement)
    cas13b_df['Target Sequence'] = cas13b_df['Guide Sequence'].apply(utils.reverse_complement)
    cas13b_df['Cas'] = 'PspCas13b'
    cas13b_df['Guide Score'] = 'Unknown'  # Placeholder for guide score, replace with actual values if available

    #### Step 2.2 - standardise Cas13d data 
    cas13d_df = utils.read_csv(cas13d_file)
    cas13d_df['Index'] = cas13d_df.index
    cas13d_df['Cas'] = 'RfxCas13d'

    
    #### Step 2.3 - combine Cas13b and Cas13d data
    common_cols = cas13b_df.columns.intersection(cas13d_df.columns).tolist()
    combined_df = pd.concat([cas13b_df[common_cols], cas13d_df[common_cols]], ignore_index=True)
    
    #### Step 2.4 - explicit order and save file
    combined_file_path = output_dir / "mers_guides.csv"
    combined_df = combined_df[['Cas', 'Index', 'Region','Guide Sequence', 'Target Sequence', 'Guide Score']]
    regions_to_keep = ["5'UTR", "3'UTR", "ORF1ab", "N", "E", "S", "M"]
    combined_df = combined_df[combined_df['Region'].isin(regions_to_keep)]

    # combined_df.to_csv(combined_file_path, index=False)

    #### Step 3 - perform blastn for guide sequences and save output
    query_fasta = utils.generate_guide_fasta(combined_df, output_dir)
    db_path = Path(reference_dir) / f"{args.accession}_subgenomic"
    # utils.blastn_query(query_fasta, db_path, output_dir, args.accession)
    # print(f"BLASTN complete. Results saved to {output_dir / (args.accession + '_blastn_output.txt')}")


    #### Step 4 - make bed files for PspCas13b and RfxCas13d guides
    blastn_output_file = output_dir / f"{args.accession}_blastn_output.txt"
    guides_file = output_dir / "mers_guides.csv"
    bed_output_dir = Path("/home/unimelb.edu.au/rbengtsson/work/library_screen/data/processed")
    utils.generate_bed_file(blastn_output_file, guides_file, bed_output_dir)




if __name__ == "__main__":
    main()