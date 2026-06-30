import subprocess
import os
import glob
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np




def make_blast_db(accession, reference_dir, concat_filename=None):
    if concat_filename is None:
        concat_filename = f"{accession}_subgenomic.fasta"

    concat_path = os.path.join(reference_dir, concat_filename)

    # Find all {accession}_*.fasta files, excluding the concat output itself
    pattern = os.path.join(reference_dir, f"{accession}_*.fasta")
    fasta_files = sorted(
        f for f in glob.glob(pattern)
        if os.path.abspath(f) != os.path.abspath(concat_path)
    )

    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found matching: {pattern}")

    print(f"Found {len(fasta_files)} FASTA files:")
    for f in fasta_files:
        print(f"  {f}")

    # Concatenate all matching FASTA files
    with open(concat_path, "w") as outfile:
        for fasta_file in fasta_files:
            with open(fasta_file, "r") as infile:
                outfile.write(infile.read())

    print(f"Concatenated to: {concat_path}")

    # Build BLAST nucleotide database
    db_path = os.path.splitext(concat_path)[0]
    subprocess.run(
        ["makeblastdb", "-in", concat_path, "-dbtype", "nucl", "-out", db_path],
        check=True,
    )
    return db_path


def blastn_query(query_fasta, db_path, output_dir, accession, pident=100, cov_perc=100):
    output_file = output_dir / f"{accession}_blastn_output.txt"
    blastn_cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", f"{query_fasta}",
        "-db", f"{db_path}",
        "-max_target_seqs", "10",
        "-perc_identity", f"{pident}",
        "-qcov_hsp_perc", f"{cov_perc}",
        "-outfmt", "6",
        "-out", f"{output_file}"]
    subprocess.run(blastn_cmd, check=True)


def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def read_csv(file_path):
    df = pd.read_csv(file_path, header=0)
    return df


def generate_guide_fasta(df, output_dir):
    records = [
        SeqRecord(Seq(seq), id=str(seq_id), description="")
        for seq_id, seq in zip(df['Cas'] + '_' + df['Index'].astype(str), df["Target Sequence"])
    ]
    SeqIO.write(records, output_dir / "query_guides.fasta", "fasta")
    return output_dir / "query_guides.fasta"


def generate_bed_file(blastn_output, mers_guides_file, output_dir):
    blastn_df = pd.read_csv(blastn_output, sep="\t", header=None,
                            names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                                   "qstart", "qend", "sstart", "send", "evalue", "bitscore"])
    
    # Map Guide sequence and target sequence to the blastn_df
    scores_df = pd.read_csv(mers_guides_file)

    blastn_df[['Cas', 'Index']] = blastn_df['qseqid'].str.split('_', n=1, expand=True)
    blastn_df['Index'] = blastn_df['Index'].astype(int)
    blastn_df['sub_genome'] = blastn_df['sseqid'].str.split('_').str[-1]
    blastn_df['sstart'] -= 1  # Convert to 0-based

    # Merge Guide Sequence and Target Sequence on Cas + Index
    seq_cols = scores_df[['Cas', 'Index', 'Guide Sequence', 'Target Sequence']]
    blastn_df = blastn_df.merge(seq_cols, on=['Cas', 'Index'], how='left')

    # Map RfxCas13d guide scores by index
    scores = scores_df[scores_df['Cas'] == 'RfxCas13d'].set_index('Index')['Guide Score']
    rfx_mask = blastn_df['Cas'] == 'RfxCas13d'
    blastn_df.loc[rfx_mask, 'Guide Score'] = blastn_df.loc[rfx_mask, 'Index'].map(scores)

    # Write BED files per Cas type
    for cas, group in blastn_df.groupby('Cas'):
        cols = ['sseqid', 'sstart', 'send', 'Index', 'sub_genome', 'Guide Sequence', 'Target Sequence']
        if cas == 'RfxCas13d':
            cols.append('Guide Score')
        group[cols].to_csv(output_dir / f"{cas}_guides.bed", sep='\t', header=False, index=False)



def visualize_guide_proximity(df, output_dir):
    """
    Visualize RfxCas13d (bed1) and PspCas13b (bed2) guide positions and proximity.

    Parameters
    ----------
    closest_df : pd.DataFrame
        Output from filter_results(). Expected columns:
        subgenomic, qstart, qend, RxfCas13d guide seq, guide score,
        PspCas13b index, sstart, send, PspCas13b guide seq, distance
    output_html : str or Path, optional
        If given, saves the figure as a self-contained HTML file.

    Returns
    -------
    fig : plotly.graph_objects.Figure
    """
    subgroups = df["subgenomic"].unique()
    n = len(subgroups)

    fig = make_subplots(
        rows=n, cols=1,
        subplot_titles=[f"{sg}" for sg in subgroups],
        shared_xaxes=False,
        vertical_spacing=0.08,
    )

    colors = {"RxfCas13d": "#1f77b4", "PspCas13b": "#ff7f0e"}
    legend_added = set()

    def assign_stagger_levels(segments, gap=0.2):
        """Assign y-level offsets so overlapping segments are staggered."""
        # segments: list of (start, end) tuples; returns list of int levels (0, 1, 2, ...)
        levels = []
        occupied = []  # list of (end, level)
        for start, end in segments:
            # find lowest free level
            used = {lvl for (occ_end, lvl) in occupied if occ_end > start}
            level = 0
            while level in used:
                level += 1
            levels.append(level)
            occupied.append((end, level))
        return levels

    for row_idx, sg in enumerate(subgroups, start=1):
        sub = df[df["subgenomic"] == sg].reset_index(drop=True)

        # Stagger RxfCas13d guides that overlap each other
        rxf_segs = list(zip(sub["qstart"], sub["qend"]))
        rxf_levels = assign_stagger_levels(rxf_segs)

        # Base y positions: PspCas13b at 1 (top), RxfCas13d base at 0 (bottom), staggered downward
        PSP_Y = 1.0
        RXF_BASE = -0.2
        STAGGER_STEP = -0.35

        for i, r in sub.iterrows():
            rxf_y = RXF_BASE + rxf_levels[i] * STAGGER_STEP

            # RxfCas13d guide
            tool = "RxfCas13d"
            show_legend = tool not in legend_added
            fig.add_trace(
                go.Scatter(
                    x=[r["qstart"], r["qend"]],
                    y=[rxf_y, rxf_y],
                    mode="lines",
                    line=dict(color=colors[tool], width=10),
                    name=tool,
                    legendgroup=tool,
                    showlegend=show_legend,
                    hovertemplate=(
                        f"<b>{tool}</b><br>"
                        f"Position: {r['qstart']}–{r['qend']}<br>"
                        f"Guide seq: {r['RxfCas13d guide seq']}<br>"
                        f"Target seq: {r['RxfCas13d target seq']}<br>"
                        f"Guide score: {r['guide score']:.4f}"
                        "<extra></extra>"
                    ),
                ),
                row=row_idx, col=1,
            )
            if show_legend:
                legend_added.add(tool)

            # PspCas13b guide
            tool = "PspCas13b"
            show_legend = tool not in legend_added
            fig.add_trace(
                go.Scatter(
                    x=[r["sstart"], r["send"]],
                    y=[PSP_Y, PSP_Y],
                    mode="lines",
                    line=dict(color=colors[tool], width=10),
                    name=tool,
                    legendgroup=tool,
                    showlegend=show_legend,
                    hovertemplate=(
                        f"<b>{tool}</b><br>"
                        f"Position: {r['sstart']}–{r['send']}<br>"
                        f"Guide seq: {r['PspCas13b guide seq']}<br>"
                        f"Target seq: {r['PspCas13b target seq']}"
                        "<extra></extra>"
                    ),
                ),
                row=row_idx, col=1,
            )
            if show_legend:
                legend_added.add(tool)

        max_level = max(rxf_levels) if rxf_levels else 0
        y_min = RXF_BASE + max_level * STAGGER_STEP - 0.4
        fig.update_yaxes(
            tickvals=[PSP_Y, RXF_BASE],
            ticktext=["PspCas13b", "RxfCas13d"],
            range=[y_min, PSP_Y + 0.4],
            row=row_idx, col=1,
        )

    fig.update_layout(
        title="Guide Position Overlap by Subgenomic Region",
        height=320 * n,
        hovermode="closest",
        legend=dict(title="Cas type", itemsizing="constant"),
        font=dict(size=12),
    )

    html_file = output_dir / "guide_overlap.html"
    fig.show()
    # fig.write_html(html_file)
    # print(f"Saved to {html_file}")