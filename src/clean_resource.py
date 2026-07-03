import pandas as pd
import sys
from pathlib import Path
from Bio import SeqIO


target_dir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(target_dir))
import src.utils as utils

# Define final schema
final_columns = ['Cas', 'Index', 'Guide Sequence', 'Target Sequence', 'Guide Score']

def get_genes_list(reference_fasta_file):
        return [record.id for record in SeqIO.parse(reference_fasta_file, "fasta")]


def clean_cas13b(file_path):
    """Standardise PspCas13b guide data."""

    columns_mapping = {
        'Spacer': 'Guide Sequence'}

    df = utils.file_to_df(file_path)
    df = df.rename(columns=columns_mapping)
    df['Index'] = df.index
    df['Target Sequence'] = df['Guide Sequence'].apply(utils.reverse_complement)
    df['Cas'] = 'PspCas13b'
    df['Guide Score'] = 'Unknown'
    df = df[final_columns]  # Reorder columns to match final schema
    return df


def clean_cas13d(file_path):
    """Standardise RfxCas13d guide data."""
    df = utils.file_to_df(file_path)
    df['Index'] = df.index
    df['Cas'] = 'RfxCas13d'
    df['Guide Score'] = df['Guide Score'].astype(float)  # Ensure guide score is float
    df = df[final_columns]  # Reorder columns to match final schema
    return df


def combine_guides(cas13b_df, cas13d_df, output_dir):
    """Merge Cas13b and Cas13d dataframes on their common columns."""
    common_cols = cas13b_df.columns.intersection(cas13d_df.columns).tolist()
    combined_df = pd.concat([cas13b_df[common_cols], cas13d_df[common_cols]], ignore_index=True)
    """save to CSV."""
    output_file = Path(output_dir) / "all_guides.csv"
    combined_df.to_csv(output_file, index=False)
    return combined_df
