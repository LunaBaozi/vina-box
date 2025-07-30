import os, argparse
import csv
import numpy as np
import pandas as pd
from matplotlib.patches import Wedge
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Add project root to path to import from scripts
project_root = Path(__file__).parent.parent.parent
sys.path.append(str(project_root))

from load_config_paths import PipelinePaths

def postprocess_vina_results(input_path, output_path):
    # Read and filter lines with at least one comma
    with open(input_path, 'r') as infile:
        lines = [line for line in infile if ',' in line]

    # Parse CSV and sort by affinity_kcal/mol
    reader = csv.DictReader(lines)
    rows = list(reader)
    rows.sort(key=lambda x: float(x['affinity_kcal/mol']))

    # Write sorted rows to output
    with open(output_path, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
        writer.writeheader()
        writer.writerows(rows)



def plot_sa_vs_affinity(sa_score_csv, vina_csv, output_plot_name, pareto_csv_path):
    # Read CSVs
    sa_df = pd.read_csv(sa_score_csv)
    vina_df = pd.read_csv(vina_csv)

    # Add the .sdf extension to sa_df['filename'] if and only if the filename is numerical
    vina_df['ligand'] = vina_df['ligand'].apply(
        lambda x: f"{x}.sdf" if str(x).isdigit() else str(x)
    )
    # Ensure both columns are string type for merging
    sa_df['filename'] = sa_df['filename'].astype(str)
    vina_df['ligand'] = vina_df['ligand'].astype(str)

    merged = pd.merge(sa_df, vina_df, left_on='filename', right_on='ligand', suffixes=('_sa', '_vina'))

    if merged.empty:
        print("Warning: No matching rows after merging. Check that 'filename' in SA score CSV matches 'ligand' in Vina CSV (with .sdf appended).")
        return

    # Prepare data
    x = merged['SA_score']
    y = merged['affinity_kcal/mol']
    # Prepare colors for SCScore (left) and NP_score (right)

    scscore = merged['SCScore']
    npscore = merged['NP_score']

    # Normalize scores for colormaps
    sc_min, sc_max = scscore.min(), scscore.max()
    np_min, np_max = npscore.min(), npscore.max()
    sc_norm = (scscore - sc_min) / (sc_max - sc_min + 1e-8)
    np_norm = (npscore - np_min) / (np_max - np_min + 1e-8)

    sc_cmap = plt.get_cmap('Oranges')  # sa score
    np_cmap = plt.get_cmap('Blues')  # np score
    sc_colors = sc_cmap(sc_norm)  #(sc_norm)
    np_colors = np_cmap(np_norm)  #(np_norm)

    # Define edge colors (default to blue, or customize as needed)
    edge_colors = merged['failed'].apply(lambda v: 'red' if pd.notna(v) and str(v).strip() else 'blue')

    # Custom scatter: draw two semicircles for each point
    fig, ax = plt.subplots(figsize=(8, 6))
    for xi, yi, sc_col, np_col, edge_col in zip(x, y, sc_colors, np_colors, edge_colors):
        # Left semicircle (SCScore)
        wedge1 = Wedge((xi, yi), 0.15, 90, 270, facecolor=sc_col, edgecolor=edge_col, linewidth=1)
        # Right semicircle (NP_score)
        wedge2 = Wedge((xi, yi), 0.15, 270, 90, facecolor=np_col, edgecolor=edge_col, linewidth=1)
        ax.add_patch(wedge1)
        ax.add_patch(wedge2)
    # Ensure circles are not stretched
    ax.set_aspect('equal', adjustable='datalim')

    ax.set_xlabel('SA_score')
    ax.set_ylabel('affinity_kcal/mol')
    ax.set_title('SA_score vs. affinity_kcal/mol')

    # Add colorbars for SCScore and NP_score, position them to the right using add_axes for consistent height
    sm_sc = ScalarMappable(cmap=sc_cmap, norm=plt.Normalize(sc_min, sc_max))
    sm_np = ScalarMappable(cmap=np_cmap, norm=plt.Normalize(np_min, np_max))

    # Shrink main plot to make space for colorbars
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])

    # Add first colorbar (SCScore)
    cbar_ax1 = fig.add_axes([0.9, box.y0, 0.02, box.height])
    cbar1 = fig.colorbar(sm_sc, cax=cbar_ax1)
    cbar1.set_label('SCScore (L)', labelpad=20, rotation=270, fontsize=10, loc='center')

    # Add second colorbar (NP_score)
    cbar_ax2 = fig.add_axes([1.05, box.y0, 0.02, box.height])
    cbar2 = fig.colorbar(sm_np, cax=cbar_ax2)
    cbar2.set_label('NP_score (R)', labelpad=20, rotation=270, fontsize=10, loc='center')

    # Add vertical lines (move labels to legend)
    ax.axvline(x=6.5, color='green', linestyle='--', linewidth=2)
    ax.axvline(x=3, color='yellow', linestyle='--', linewidth=2)
    ax.axvline(x=2, color='orange', linestyle='--', linewidth=2)

    # Add horizontal lines (move labels to legend)
    ax.axhline(y=-6, color='pink', linestyle=':', linewidth=2)
    ax.axhline(y=-8, color='purple', linestyle=':', linewidth=2)
    ax.axhline(y=-9, color='blue', linestyle=':', linewidth=2)

    # Find indices of points to label
    idx_affinity = merged['affinity_kcal/mol'].astype(float).idxmin()
    idx_sa = merged['SA_score'].astype(float).idxmin()
    idx_sc = merged['SCScore'].astype(float).idxmin()
    idx_np = merged['NP_score'].astype(float).idxmin()

    label_indices = set([idx_affinity, idx_sa, idx_sc, idx_np])

    for idx in label_indices:
        xi = merged.loc[idx, 'SA_score']
        yi = merged.loc[idx, 'affinity_kcal/mol']
        ligand_label = merged.loc[idx, 'ligand']
        ax.annotate(
            ligand_label,
            (xi, yi),
            textcoords="offset points",
            xytext=(0, 15),  # Place label below the point
            ha='center',
            va='top',
            fontsize=8,
            color='black',
            bbox=dict(boxstyle="round,pad=0.2", fc="white", alpha=0.7, lw=0)
        )

    

    # --- Pareto frontier (minimize both SA_score and affinity_kcal/mol) ---
    # Sort by SA_score, then iterate to keep only points with lowest y so far
    pareto_points = merged.sort_values(['SA_score', 'affinity_kcal/mol'])
    pareto_front = []
    pareto_indices = []
    min_y = float('inf')
    for idx, row in pareto_points.iterrows():
        if row['affinity_kcal/mol'] < min_y:
            pareto_front.append((row['SA_score'], row['affinity_kcal/mol']))
            pareto_indices.append(idx)
            min_y = row['affinity_kcal/mol']
    if pareto_front:
        pareto_front = np.array(pareto_front)
        ax.plot(pareto_front[:, 0], pareto_front[:, 1], color='green', linewidth=2, marker='o', markersize=4, label='Pareto frontier')
        # Save Pareto front rows to CSV
        pareto_df = merged.loc[pareto_indices]
        os.makedirs(os.path.dirname(pareto_csv_path), exist_ok=True)
        pareto_df.to_csv(pareto_csv_path, index=False)
        # Add to legend
        handles, labels = ax.get_legend_handles_labels()
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='red', markersize=10, label='Failed'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='blue', markersize=10, label='Passed'),
            Line2D([0], [0], color='green', linewidth=2, marker='o', markersize=6, label='Pareto frontier')
        ]
        labels.append('Pareto frontier')
        ax.legend(handles=legend_elements, loc='upper left', title='Lipinski rules')
    else:  # Add a legend for the border colors of the sample points
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='red', markersize=10, label='Failed'),
            Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markeredgecolor='blue', markersize=10, label='Passed')
        ]
        ax.legend(handles=legend_elements, loc='upper left', title='Lipinski rules')

    plt.savefig(output_plot_name, bbox_inches='tight')
    plt.show()
    




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Post-process Vina docking results and create analysis plots.")
    parser.add_argument('--num_gen', type=int, required=False, default=0, help='Desired number of generated molecules (int, positive)')
    parser.add_argument('--epoch', type=int, required=False, default=0, help='Epoch number the model will use to generate molecules (int, 0-99)')
    parser.add_argument('--known_binding_site', type=str, required=False, default='0', help='Allow model to use binding site information (True, False)')
    parser.add_argument('--aurora', type=str, required=False, default='B', help='Aurora kinase type (str, A, B)')
    parser.add_argument('--pdbid', type=str, required=False, default='4af3', help='PDB ID (str)')
    parser.add_argument('--experiment', type=str, required=False, default='default', help='Experiment name (str)')
    args = parser.parse_args()

    # Initialize paths
    paths = PipelinePaths()

    epoch = args.epoch
    num_gen = args.num_gen
    known_binding_site = args.known_binding_site
    aurora = args.aurora
    pdbid = args.pdbid.lower()
    experiment = args.experiment

    # Use PipelinePaths to get consistent directory structure
    vina_base_dir = paths.vina_box_docking_path(pdbid, experiment, epoch, num_gen, known_binding_site, pdbid)
    
    # Input files
    vina_csv = os.path.join(vina_base_dir, "vina_results.csv")
    synth_csv = paths.hope_box_results_path(experiment, epoch, num_gen, known_binding_site, pdbid, 'merged_scores.csv')
    
    # Output files
    vina_postprocess = os.path.join(vina_base_dir, "vina_results_postprocessed.csv")
    output_plot_name = os.path.join(vina_base_dir, "sa_vs_affinity_plot.png")
    pareto_csv_path = os.path.join(vina_base_dir, "pareto_front.csv")

    # Create output directories
    os.makedirs(os.path.dirname(vina_postprocess), exist_ok=True)
    os.makedirs(os.path.dirname(output_plot_name), exist_ok=True)

    # Check if input files exist
    if not os.path.exists(vina_csv):
        print(f"Error: Vina results file not found: {vina_csv}")
        exit(1)
    
    if not os.path.exists(synth_csv):
        print(f"Error: Synthesizability scores file not found: {synth_csv}")
        exit(1)

    print(f"Processing Vina results from: {vina_csv}")
    print(f"Using synthesizability scores from: {synth_csv}")

    postprocess_vina_results(
        vina_csv,
        vina_postprocess
    )
    
    plot_sa_vs_affinity(
        synth_csv,
        vina_postprocess,
        output_plot_name,
        pareto_csv_path
    )

    print(f"Results saved:")
    print(f"  Postprocessed Vina results: {vina_postprocess}")
    print(f"  SA vs Affinity plot: {output_plot_name}")
    print(f"  Pareto front: {pareto_csv_path}")
    print("Analysis complete!")

    