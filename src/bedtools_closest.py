import os
from Bio import SeqIO
from pathlib import Path
import sys
import argparse
import subprocess
import pybedtools
import utils as utils
import visualisation as visualisation



def parse_args():
    parser = argparse.ArgumentParser(description="Find closest features between BED files")

    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument("-a", 
                          type=str, required=True,
                          help="PspCas13b BED file", )
    required.add_argument("-b", 
                          type=str, required=True,
                          help="RfxCas13d BED file")
    required.add_argument("-o", 
                          type=str, required=True,
                          help="Output directory")
    
    
    # Optional arguments
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-f", 
                          type=str, default=None,
                          help="Output file name")
    optional.add_argument("-w", 
                          type=int, default=None,
                          help="Filter by window size in bp (50bp default)")
    return parser.parse_args()
    

def run_bedtools_srt(query_file, subject_file):
    query = pybedtools.BedTool(query_file)
    subject = pybedtools.BedTool(subject_file)

    query_sorted = query.sort()
    subject_sorted = subject.sort()

    return query_sorted, subject_sorted


def filter_guide_score(query_file, threshold=0.8):
    query = pybedtools.BedTool(query_file)
    return query.filter(lambda x: float(x[6]) > threshold)  # Assuming guide score is the 8th column


def run_bedtools_closest(query_sorted, subject_sorted):
    closest = query_sorted.closest(subject_sorted, t="all", d=True)
    return closest



def filter_results(closest_features):

    # Convert to DataFrame for easier manipulation
    df = closest_features.to_dataframe(names=["gene", "qstart", "qend", "PspCas13b index", 
                                              "PspCas13b guide seq", "PspCas13b target seq",
                                              "closest_chrom", "sstart", "send",
                                              "RfxCas13d index", "RfxCas13d guide seq", "RfxCas13d target seq", 
                                              "guide score", "distance"])
    

    # Drop columns that are not needed for the final output
    df = df.drop(columns=["closest_chrom"])

    # Change column order to match the desired output
    df = df[["gene", "PspCas13b index", "qstart", "qend", "PspCas13b guide seq", "PspCas13b target seq",
             "RfxCas13d index", "sstart", "send",
             "RfxCas13d guide seq", "RfxCas13d target seq", "guide score", "distance"]]  
    
    # select top 10 for each group based on distance and score
    df_filt = (df.sort_values(by=["distance", "guide score"], ascending=[True, False])
                 .groupby("gene")
                 .head(10)
                 .reset_index(drop=True))

    return df_filt



def main():
    args = parse_args()
    bed_file1 = Path(args.a) # PspCas13b BED file
    bed_file2 = Path(args.b) # RfxCas13d BED file
    output_dir = Path(args.o)

    ### Sort the BED files
    sorted_bed1, sorted_bed2 = run_bedtools_srt(bed_file1, bed_file2)


    ### Filter the first BED file based on guide score
    filt_srt_bed2 = filter_guide_score(sorted_bed2, threshold=0.8)

    if args.w is not None:
        window_size = args.w
        nearby = sorted_bed1.window(filt_srt_bed2, w=window_size)
        filtered_results = filter_results(nearby)
        
        ### Save the filtered results to a CSV file
        bedfile_dir = Path(os.path.dirname(bed_file1))
        output_file = bedfile_dir / f"closest_guides_{args.w}bp.csv"
        # print(f"Saving filtered results to {output_file}")
        filtered_results.to_csv(output_file, index=False)

    
    else:
        ### Get closest coordinates
        closest_features = run_bedtools_closest(sorted_bed1, filt_srt_bed2)

        ### Sort the closest features by distance (last column and score)
        filtered_results = filter_results(closest_features)
    
        ### Save the filtered results to a CSV file
        bedfile_dir = Path(os.path.dirname(bed_file1))
        output_file = bedfile_dir / "closest_guides.csv"
        filtered_results.to_csv(output_file, index=False)


    # ### Visualize the results
    # outfile_name = args.f
    # fig = visualisation.visualize_guide_proximity(filtered_results, output_dir, bed_file1, outfile_name)


if __name__ == "__main__":
    main()