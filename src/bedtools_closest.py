from Bio import SeqIO
from pathlib import Path
import sys
import argparse
import subprocess
import pybedtools
import utils



def parse_args():
    parser = argparse.ArgumentParser(description="Find closest features between BED files")
    parser.add_argument("-a", help="Query BED file", type=str, required=True)
    parser.add_argument("-b", help="Subject BED file", type=str, required=True)
    parser.add_argument("-o", help="Output directory", type=str, required=True)
    return parser.parse_args()
    

def run_bedtools_srt(query_file, subject_file):
    query = pybedtools.BedTool(query_file)
    subject = pybedtools.BedTool(subject_file)

    query_sorted = query.sort()
    subject_sorted = subject.sort()

    return query_sorted, subject_sorted


def filter_guide_score(query_file, threshold=0.8):
    query = pybedtools.BedTool(query_file)
    return query.filter(lambda x: float(x[7]) > threshold)  # Assuming guide score is the 8th column


def run_bedtools_closest(query_sorted, subject_sorted):
    closest = query_sorted.closest(subject_sorted, d=True)
    return closest



def filter_results(closest_features):

    # Convert to DataFrame for easier manipulation
    df = closest_features.to_dataframe(names=["chrom", "qstart", "qend", "RxfCas13d index", "subgenomic", 
                                              "RxfCas13d guide seq", "RxfCas13d target seq", "guide score",
                                              "closest_chrom", "sstart", "send",
                                              "PspCas13b index", "subgenomic 2",
                                              "PspCas13b guide seq", "PspCas13b target seq", "distance"])
    
    # Drop columns that are not needed for the final output
    df = df.drop(columns=["chrom", "closest_chrom", "subgenomic 2"])

    # Change column order to match the desired output
    df = df[["subgenomic", "qstart", "qend", "RxfCas13d index", 
             "RxfCas13d guide seq", "RxfCas13d target seq", "guide score", 
             "PspCas13b index", "sstart", "send",
             "PspCas13b guide seq", "PspCas13b target seq", "distance"]]
    
    
    # select top 10 for each group based on distance and score
    df_filt = (df.sort_values(by=["distance", "guide score"], ascending=[True, False])
                 .groupby("subgenomic")
                 .head(10)
                 .reset_index(drop=True))

    return df_filt



def main():
    args = parse_args()
    bed_file1 = Path(args.a)
    bed_file2 = Path(args.b)
    output_dir = Path(args.o)

    ### Sort the BED files
    sorted_bed1, sorted_bed2 = run_bedtools_srt(bed_file1, bed_file2)

    ### Filter the first BED file based on guide score
    filt_srt_bed1 = filter_guide_score(sorted_bed1, threshold=0.8)

    ### Get closest coordinates
    closest_features = run_bedtools_closest(filt_srt_bed1, sorted_bed2)

    ### Sort the closest features by distance (last column and score)
    filtered_results = filter_results(closest_features)

    ### Save the filtered results to a CSV file
    output_file = output_dir / "closest_guides.csv"
    # filtered_results.to_csv(output_file, index=False)

    ### Visualize the results
    fig = utils.visualize_guide_proximity(filtered_results, output_dir)
    # fig.show()


if __name__ == "__main__":
    main()