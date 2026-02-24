#!/usr/bin/env python3

import pysam
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import re
from matplotlib.backends.backend_pdf import PdfPages
import argparse
from collections import defaultdict

# -----------------------
# ARGUMENTS
# -----------------------
parser = argparse.ArgumentParser(description="Barcode QC report generator")
parser.add_argument("--bam_folder", type=str, required=True, help="Path to folder containing BAM files")
parser.add_argument("--barcode_folder", type=str, required=True, help="Path to folder containing expected barcodes CSVs")
parser.add_argument("--outputfile", type=str, required=True, help="filenename to save PDF report")
parser.add_argument("--numcells", type=int, required=False, help="number of cells sorted")

args = parser.parse_args()

bam_folder = Path(args.bam_folder)
barcode_folder = Path(args.barcode_folder)
outputfile = Path(args.outputfile)
numcells = args.numcells

tags = ["XD", "XE", "XF"]  # barcode tags
umi_tag = "XM"              # UMI tag
xc_tag = "XC"               # 3-barcode combination tag

# Load valid barcodes into a set
valid_combos = barcode_folder / "all_combinations.txt"
with open(valid_combos) as f:
    valid_xc = set(line.strip() for line in f if line.strip())

rows = list("ABCDEFGH")
cols = list(range(1, 13))

# -----------------------
# FUNCTIONS
# -----------------------
def read_expected_barcodes(file_path):
    df = pd.read_csv(file_path, header=0, dtype=str)
    df = df.rename(columns={"WellPosition": "well", "Barcode": "barcode"})
    df['barcode'] = df['barcode'].str.strip()
    df['well'] = df['well'].str.strip()
    return df.sort_values("barcode")

def extract_all_counts(bam_file, tags=["XD","XE","XF"], xc_tag="XC", umi_tag="XM"):
    observed_counts = {tag: Counter() for tag in tags}
    umi_counts = Counter()
    combo_counts = Counter()
    umi_counts_per_combo = Counter()
    umi_seq_counts = Counter()
    wrong_umi_count = 0
    missing_barcode_count = 0
    total_reads = 0

    with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as bam:
        for read in bam:
            total_reads += 1

            # Barcode counts
            for tag in tags:
                if read.has_tag(tag):
                    observed_counts[tag][read.get_tag(tag)] += 1

            # Check UMI and barcode combination
            xc_val = read.get_tag(xc_tag) if read.has_tag(xc_tag) else None
            xm_val = read.get_tag(umi_tag) if read.has_tag(umi_tag) else None

            if xm_val is None or len(xm_val) != 10:
                wrong_umi_count += 1

            if xc_val is None or len(xc_val) < 18 or xc_val not in valid_xc:
                missing_barcode_count += 1

            # Only count valid UMI combinations
            if xc_val is not None and xm_val is not None and xc_val in valid_xc:
                key = f"{xc_val}_{xm_val}"
                umi_counts[key] += 1

            # only count valid UMI_seq combinations = UMI tag + combination bc1, bc2 and bc3 + sequence
            if xc_val is not None and xm_val is not None and xc_val in valid_xc:
                seq_val = read.query_sequence  # or read.get_forward_sequence() if soft-clipped
                if seq_val:
                    key_seq = f"{xc_val}_{xm_val}_{seq_val}"
                    umi_seq_counts[key_seq] += 1
            
            # Count reads per combination (unique XC)
            if xc_val is not None and xc_val in valid_xc:
                combo_counts[xc_val] += 1

    # Estimation of cell number
    # number of unique XC that have at least 20 reads
    for key, count in umi_seq_counts.items():
        xc_val, xm_val, seq_val = key.split("_", 2)  # split only on the first "_"
        umi_counts_per_combo[xc_val] += 1   # count number of distinct UMIseqs
    est_cell = sum(1 for xc, count in umi_counts_per_combo.items() if count >= 20)

    # Compute UMI collision metrics
    seqs_per_umi = defaultdict(set)
    for key in umi_seq_counts.keys():
        xc_val, xm_val, seq_val = key.split("_", 2)
        seqs_per_umi[(xc_val, xm_val)].add(seq_val)

    num_umis_with_multiple_seqs = sum(len(s) > 1 for s in seqs_per_umi.values())
    mean_num_seqs_per_umi = np.mean([len(s) for s in seqs_per_umi.values()]) if seqs_per_umi else 0

    # Convert to DataFrames
    observed_dfs = {tag: pd.DataFrame(counter.items(), columns=["barcode","count"]).sort_values("barcode")
                    for tag, counter in observed_counts.items()}

    #return observed_dfs, combo_counts, est_cell, umi_counts, wrong_umi_count, missing_barcode_count, total_reads
    return observed_dfs, combo_counts, est_cell, umi_counts, wrong_umi_count, missing_barcode_count, total_reads, umi_seq_counts, num_umis_with_multiple_seqs, mean_num_seqs_per_umi

def join_with_expected(observed_df, expected_df):
    merged = pd.merge(observed_df, expected_df, on="barcode", how="outer", indicator=True)
    counts = merged["_merge"].value_counts()
    num_correct_barcodes = counts.get("both", 0)          # number of barcodes found in both
    num_wrong_barcodes = counts.get("left_only", 0)     # in observed_df but not expected_df
    num_missing_barcodes = counts.get("right_only", 0)   # in expected_df but not observed_df
    return merged, num_correct_barcodes, num_wrong_barcodes, num_missing_barcodes

def parse_well_position(well):
    match = re.match(r"([A-H])(\d{1,2})", well)
    if match:
        return match.group(1), int(match.group(2))
    return None, None

def plot_pie_chart(counts, lib_name, ax):
    sizes = [
        sum(v for v in counts.values() if v == 1),
        sum(v for v in counts.values() if 1 < v <= 25),
        sum(v for v in counts.values() if 25 < v <= 50),
        sum(v for v in counts.values() if 50 < v <= 100),
        sum(v for v in counts.values() if v > 100),
    ]
    labels = ["Count = 1", "1 < Count ≤ 25", "25 < Count ≤ 50", "50 < Count ≤ 100", "Count > 100"]
    colors = ["green", "lightblue", "gold", "pink", "lightcoral"]

    ax.clear()
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors, startangle=140)
    ax.set_title(f"Molecule count distribution: {lib_name}")

def plot_heatmap(df, lib_name, bc_name, ax, tag_map):
    df[['row', 'col']] = df['well'].apply(lambda x: pd.Series(parse_well_position(x)))
    heatmap_data = df.pivot(index="row", columns="col", values="count").reindex(index=rows, columns=cols)

    ax.clear()
    sns.heatmap(
        heatmap_data,
        cmap="viridis",
        cbar=True,
        vmin=heatmap_data.min().min(),
        vmax=heatmap_data.max().max(),
        ax=ax,
        square=True,  # force equal aspect ratio
        cbar_kws={"shrink": 0.6}  # make the colorbar 60% of the height
    )

    label = tag_map.get(bc_name, bc_name)
    ax.set_title(f"{lib_name} {label} ({bc_name})")
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")

def make_stats_table(stats_df, pdf, n_panels=4):
    """
    Draw summary table on a page with same dimensions as library pages.
    n_panels: number of panels on library pages (to match width)
    """
    stats_df_T = stats_df.T

    fig, ax = plt.subplots(figsize=(5 * n_panels, 5))  # fixed width/height
    ax.axis("off")

    cell_text = stats_df_T.astype(str).values.tolist()

    table = ax.table(
        cellText=cell_text,
        rowLabels=stats_df_T.index.tolist(),
        colLabels=None,
        loc="center"
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.2)  # vertical scaling

    ax.set_title("Library Statistics", fontsize=14, weight="bold")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def make_library_page(lib_name, umi_counts, merged_dfs, pdf, tag_map):
    n_panels = len(merged_dfs) + 1
    fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 5))

    if isinstance(axes, np.ndarray):
        axes = axes.flatten().tolist()
    else:
        axes = [axes]

    # Pie chart
    #plot_pie_chart(umi_counts, lib_name, axes[0])
    plot_pie_chart(umi_counts, lib_name, axes[0])

    # Heatmaps
    for i, (tag, df) in enumerate(merged_dfs.items(), start=1):
        plot_heatmap(df, lib_name, tag, axes[i], tag_map)

    # Hide unused axes
    for j in range(len(merged_dfs) + 1, len(axes)):
        axes[j].axis("off")

    fig.suptitle(f"Library: {lib_name}", fontsize=16, weight="bold")
    plt.tight_layout(rect=[0, 0, 1, 0.97])
    pdf.savefig(fig)
    plt.close(fig)

# -----------------------
# MAIN LOOP
# -----------------------

# First pass: collect stats only
bam_files = list(bam_folder.glob("*_tagged.bam"))
#bam_files = list(bam_folder.glob("*_tagging.bam"))
all_stats = []
library_data = []

for bam_file in bam_files:
    lib_base = bam_file.stem.replace("_tagged", "")
    print(f"Processing {lib_base}")
    #observed_dfs, combo_counts, est_cell, umi_counts, wrong_umi, missing_barcode_count, total_reads = extract_all_counts(bam_file)
    (observed_dfs, combo_counts, est_cell, umi_counts, wrong_umi, missing_barcode_count, total_reads,umi_seq_counts, num_umis_with_multiple_seqs, mean_num_seqs_per_umi) = extract_all_counts(bam_file)
    
    # Saving the umiseq stats
    df_umiseq_counts = pd.DataFrame(
    [(k, v) for k, v in umi_seq_counts.items()],
    columns=["key_seq", "count"]
    )
    out_file = f"{lib_base}_umi_seq_counts.tsv"
    df_umiseq_counts.to_csv(out_file, sep="\t", index=False)

    merged_dfs = {}
    wrong_barcode_df = {}
    missing_in_observed_df = {}
    for i, tag in enumerate(tags, start=1):
        expected_file = barcode_folder / f"expected_barcodes_{i}.csv"
        expected_df = read_expected_barcodes(expected_file)
        merged_df, num_correct_barcodes, num_wrong_barcodes, num_missing_barcodes = join_with_expected(observed_dfs[tag], expected_df)
        merged_dfs[tag] = merged_df
        wrong_barcode_df[tag] = num_wrong_barcodes
        missing_in_observed_df[tag] = num_missing_barcodes

    #umi_vals = list(umi_counts.values())
    umi_seq_vals = list(umi_seq_counts.values())

    if numcells is not None:
        reads_per_cell = round(total_reads / numcells)
    else:
        reads_per_cell = "NA" 

    stats_row = {
        "Library": lib_base,
        "Total Reads": total_reads,
        "# Barcode Combinations": len(combo_counts),
        "Combinations with at least 20 UMIs": est_cell, 
        #"# wrong Barcode 1": wrong_barcode_df["XD"],
        #"# wrong Barcode 2": wrong_barcode_df["XE"],
        #"# wrong Barcode 3": wrong_barcode_df["XF"],
        #"# UMItag + barcode combination": len(umi_counts),
        "# UMIs": len(umi_seq_counts),
        #"# UMIs / combination": len(umi_seq_vals)/len(combo_counts),
        "# UMIs / combination (min 20UMIs)": len(umi_seq_vals)/len(combo_counts),
        "# Reads missing or wrong length UMI (%)": f"{wrong_umi} ({wrong_umi/total_reads:.2%})",
        "Saturation (1 - (#UMI / #total))": 1-(len(umi_seq_counts)/total_reads),
        "# reads / cell": round(reads_per_cell),
        #"# Reads missing or with invalid barcode combo (%)": f"{missing_barcode_count} ({missing_barcode_count/total_reads:.2%})",
        #"Mean UMI duplication": round(np.mean(umi_vals),2),
        #"Median UMI duplication": round(np.median(umi_vals),2),
        #"Min number of UMI": np.min(umi_vals),
        #"Max Number of UMI": np.max(umi_vals),
        #"Std_dev of UMI duplication": round(np.std(umi_vals),2),
        #"UMIs with ≥2 distinct sequences": num_umis_with_multiple_seqs,
        "Mean # sequences per UMI": round(mean_num_seqs_per_umi, 2),
        "Mean UMI duplication": round(np.mean(umi_seq_vals),2),
        "Median UMI duplication": round(np.median(umi_seq_vals),2),
        "Min number of UMI": np.min(umi_seq_vals),
        "Max Number of UMI": np.max(umi_seq_vals),
        "Std_dev of UMI duplication": round(np.std(umi_seq_vals),2)
    }
    all_stats.append(stats_row)
    library_data.append((lib_base, umi_seq_counts, merged_dfs))

# Second pass: write PDF
pdf_out = outputfile
tag_map = {"XD": "Barcode 1", "XE": "Barcode 2", "XF": "Barcode 3"}
print(f"Creating Report")
with PdfPages(pdf_out) as pdf:
    # Summary page first
    stats_df = pd.DataFrame(all_stats)
    print("Creating Summary Table")
    make_stats_table(stats_df, pdf)

    # Library pages
    print("Creating plots for each library")
    for lib_base, umi_seq_counts, merged_dfs in library_data:
        make_library_page(lib_base, umi_seq_counts, merged_dfs, pdf, tag_map)


print(f"Report saved to {pdf_out}")
print("All done!")

