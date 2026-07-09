"""
This script provides functions for visualizing guide RNA positions for Cas13b and Cas13d systems using Plotly. 
It generates interactive plots showing the overlap of guide positions across different subgenomic regions.
"""

import os, sys
from zipfile import Path
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import entrez as entrez
from pathlib import Path
import numpy as np


def get_gene_length(gene_coords_dict):
    for gene, coordinates in gene_coords_dict.items():
        start, end = coordinates["start"], coordinates["end"]
        gene_lengths = int(end) - int(start)
        gene_coords_dict[gene] = {"start": start, "end": end, "length": gene_lengths}
    return gene_coords_dict



def visualize_guide_proximity(df, output_dir, bed_file_path, 
                              outfile_name="closest_features.html"):

    # Get gene coordinates from Entrez
    acession = Path(bed_file_path).parts[-2]

    gbk_dir_path = os.path.dirname(bed_file_path)
    gbk_file_path = os.path.join(gbk_dir_path, f"{acession}.gbk")
    gene_coords = entrez.get_feature_locations(gbk_file_path)

    gene_coords = get_gene_length(gene_coords)


    # Create a subplot for each unique gene
    genes = df["gene"].unique()
    genes = [g for g in genes if g != "orf1ab_1"]
    n = len(genes) 

    def assign_tracks(df, start_col='qstart', end_col='qend'):
        """
        Assign a track (y-level) to each interval so overlapping guides
        are stacked, non-overlapping guides share a track.
        """
        df = df.sort_values(start_col).reset_index(drop=True)
        track_ends = []  # tracks[i] = end position of last interval placed on that track
        
        tracks = []
        for _, row in df.iterrows():
            placed = False
            for t, end in enumerate(track_ends):
                if row[start_col] >= end:  # no overlap with last item on this track
                    track_ends[t] = row[end_col]
                    tracks.append(t)
                    placed = True
                    break
            if not placed:
                track_ends.append(row[end_col])
                tracks.append(len(track_ends) - 1)

        df['track'] = tracks
        return df


    # Pre-compute track counts so row heights scale with content
    track_counts = []
    for gene in genes:
        gene_df = df[df["gene"] == gene].reset_index(drop=True)
        pspcas13b = gene_df.iloc[:,0:6].groupby('PspCas13b index').first().reset_index()
        rfxcas13d = (gene_df.iloc[:, np.r_[0, 6:12]]
                    .groupby('RfxCas13d index').first().reset_index())
        psp_n = assign_tracks(pspcas13b.copy(), start_col='qstart', end_col='qend')['track'].nunique()
        rfx_n = assign_tracks(rfxcas13d.copy(), start_col='sstart', end_col='send')['track'].nunique()
        track_counts.append(psp_n + rfx_n + 1)  # +1 for separator gap

    total_tracks = sum(track_counts)
    row_heights = [tc / total_tracks for tc in track_counts]


    ##### Generate plots
    Y_PAD = 2

    colors = {"RxfCas13d": "#1f77b4", "PspCas13b": "#ff7f0e"}
    legend_added = set()

    fig = make_subplots(
        rows=n, cols=1,
        subplot_titles=[f"{gene}" for gene in genes],
        shared_xaxes=False,
        vertical_spacing=0.05,
        row_heights=row_heights,
    )

    psp_track_ns = {}
    rfx_track_ns = {}

    for row_idx, gene in enumerate(genes, start=1):

        gene_df = df[df["gene"] == gene].reset_index(drop=True)
        pspcas13b = gene_df.iloc[:,0:6].groupby('PspCas13b index').first().reset_index()
        rfxcas13d = (gene_df.iloc[:, np.r_[0, 6:12]]
                    .groupby('RfxCas13d index').first().reset_index())

        # PspCas13b — positive integer tracks (top)
        psp_sub = assign_tracks(pspcas13b, start_col='qstart', end_col='qend')
        psp_track_ns[row_idx] = psp_sub['track'].nunique()

        cas_type = "PspCas13b"
        show_legend = cas_type not in legend_added
        legend_added.add(cas_type)
        for _, row in psp_sub.iterrows():
            fig.add_trace(go.Scatter(
                x=[row['qstart'], row['qend']],
                y=[row['track'], row['track']],
                mode='lines',
                line=dict(width=10, color=colors["PspCas13b"]),
                name=cas_type,
                legendgroup=cas_type,
                showlegend=show_legend,
                hovertemplate=(
                    f"<b>PspCas13b</b><br>"
                    f"Position: {row['qstart']} - {row['qend']}<br>"
                    f"Guide seq: {row['PspCas13b guide seq']}<br>"
                    "<extra></extra>"
                ),
            ),
            row=row_idx, col=1)
            show_legend = False  # only first trace shows in legend

        # RfxCas13d — negative integer tracks (bottom)
        rfx_sub = assign_tracks(rfxcas13d, start_col='sstart', end_col='send')
        rfx_track_ns[row_idx] = rfx_sub['track'].nunique()
        cas_type = "RxfCas13d"
        show_legend = cas_type not in legend_added
        legend_added.add(cas_type)
        for _, row in rfx_sub.iterrows():
            fig.add_trace(go.Scatter(
                x=[row['sstart'], row['send']],
                y=[-(row['track'] + 1), -(row['track'] + 1)],
                mode='lines',
                line=dict(width=10, color=colors["RxfCas13d"]),
                name=cas_type,
                legendgroup=cas_type,
                showlegend=show_legend,
                hovertemplate=(
                    f"<b>RfxCas13d</b><br>"
                    f"Guide idx: {row['RfxCas13d index']}<br>"
                    f"Position: {row['sstart']} - {row['send']}<br>"
                    f"Guide seq: {row['RfxCas13d guide seq']}<br>"
                    f"Target seq: {row['RfxCas13d target seq']}<br>"
                    f"Guide score: <b>{row['guide score']:.3f}</b><br>"
                    "<extra></extra>"
                ),
            ),
            row=row_idx, col=1)
            show_legend = False  # only first trace shows in legend

        ymax = psp_track_ns[row_idx] - 0.5 + Y_PAD
        ymin = -rfx_track_ns[row_idx] - 0.5 - Y_PAD
        fig.update_yaxes(range=[ymin, ymax], row=row_idx, col=1)

        x_end = gene_coords[gene]["length"]
        fig.update_xaxes(range=[0, x_end], row=row_idx, col=1)


    fig.update_xaxes(title_text='Position (nt)', row=n, col=1)
    fig.update_yaxes(showticklabels=False, title_text=None)

    fig.update_layout(
        title="Guide position overlap by subgenomic region",
        autosize=True,
        height=320 * n,
        # width=900,
        margin=dict(l=30, r=30, t=60, b=20),
        legend=dict(title="Cas type", itemsizing="constant"),
    )

    html_file = output_dir / outfile_name
    # fig.show()
    fig.write_html(html_file)
    print(f"Saved to {html_file}")