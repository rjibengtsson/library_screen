import os, sys
from Bio import SeqIO
from pathlib import Path
import argparse
from numpy import indices
import pandas as pd
import typing



def parse_args():
    parser = argparse.ArgumentParser(description="Get RfxCas13d guide sequences from csv file")
    
    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '-csv', '--csv_file',
        type=str, required=True,
        help='Input CSV file containing RfxCas13d and PspCas13b bedtool closest results'
    )
    required.add_argument(
        '-in', '--input',
        type=str, required=True,
        help='Text file containing list of RfxCas13d indices to extract guide sequences for'
    )
    required.add_argument(
        '-o', '--output_file',
        type=str, required=True,
        help='Output file to write the extracted guide sequences'
    )

    return parser.parse_args()



def get_guide_sequences(csv_file: Path, input_file: Path, output_file: Path) -> typing.Dict[str, str]:
    
    # Read the list of indices from the input text file
    with open(input_file, 'r') as f:
        indices = [int(line.strip()) for line in f.readlines()]

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file, header=0)
    
    # Filter the DataFrame to include only the specified indices
    filtered_df = df[df['RfxCas13d index'].isin(indices)]

    # Wrtie output to a csv file
    result_df = filtered_df[['RfxCas13d index', 'RfxCas13d guide seq', 
                             'PspCas13b guide seq', 'gene', 'guide score']]

    columns_mapping = {
    'RfxCas13d index': 'ID',
    'RfxCas13d guide seq': 'Spacer',
    'PspCas13b guide seq': 'Closest PspCas13b spacer',
    'gene': 'Gene',
    'guide score': 'TIGER predicted guide score'
    }

    result_df.rename(columns=columns_mapping)
    return result_df.to_csv(output_file, index=False)




def main():
    args = parse_args()
    csv_file = Path(args.csv_file)
    input_file = Path(args.input)
    output_file = Path(args.output_file)
    # Get the guide sequences
    get_guide_sequences(csv_file, input_file, output_file)


if __name__ == "__main__":
    main()