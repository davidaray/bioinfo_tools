import argparse
import gzip
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

# Define TE classes and colors
TE_CLASSES = ["LINE", "SINE", "LTR", "DIRS", "DNA", "RC", "Unknown", "Satellite", "Other", "NonLTR", "Unmasked"]
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22",
          "#17becf", "#aec7e8", "#000000"]

def parse_summary_file(summary_file):
    """Extract genome size from the summary.gz file."""
    with gzip.open(summary_file, 'rt') as f:
        for line in f:
            if line.startswith("Total Length:"):
                return int(line.split(":")[1].strip().split()[0])

def parse_repeatmasker_out(rm_file, genome_size, bin_size, classes):
    """Parse the RepeatMasker .out.gz file and calculate TE class proportions."""
    with gzip.open(rm_file, 'rt') as f:
        lines = f.readlines()[3:]  # Skip header
    data = []
    for line in lines:
        fields = line.split()
        if len(fields) < 11:
            continue
        divergence = float(fields[1])
        length = abs(int(fields[6]) - int(fields[5]))
        te_class = fields[10].split('/')[0]
        if te_class == "Simple_repeat":
            te_class = "Satellite"
        elif te_class not in TE_CLASSES:
            te_class = "Other"
        if te_class in classes:
            data.append((divergence, length, te_class))

    df = pd.DataFrame(data, columns=["Divergence", "Length", "Class"])
    df["Bin"] = (df["Divergence"] // bin_size).astype(int) * bin_size
    df["Proportion"] = df["Length"] / genome_size
    return df.groupby(["Bin", "Class"]).sum().reset_index()

def create_stacked_bar_plot(df, classes, colors, max_divergence, spacing, output_prefix):
    """Generate a stacked bar plot."""
    pivot = df.pivot(index="Bin", columns="Class", values="Proportion").fillna(0)
    bins = pivot.index
    
    # Output the pivoted data to a .tsv file
    output_file = f"{output_prefix}_stacked_bar_data.tsv"
    pivot.to_csv(output_file, sep='\t', index_label='Bin')
    print(f"Data has been written to {output_file}")
    
    plt.figure(figsize=(10, 6))
    bottom = 0
    for i, cls in enumerate(classes):
        if cls in pivot.columns:
            plt.bar(
                bins, 
                pivot[cls], 
                bottom=bottom, 
                label=cls, 
                color=colors[TE_CLASSES.index(cls)],
                edgecolor="none"
            )
            bottom += pivot[cls]
    plt.title("Proportion by Genetic Divergence/Class", fontsize=14)
    plt.xlabel("Genetic Divergence (%)")
    plt.ylabel("Proportion of Genome")
    # Adjusted x-axis limits for reduced offset and full visibility of the first bar
    plt.xlim(-spacing * 0.1, max_divergence + spacing / 2)
    plt.legend(loc="upper right")
    plt.savefig(f"{output_prefix}_stackedbar.png", bbox_inches="tight")
    plt.close()

def create_line_plot(df, classes, colors, max_divergence, spacing, output_prefix):
    """Generate a line plot."""
    pivot = df.pivot(index="Bin", columns="Class", values="Proportion").fillna(0)
    bins = pivot.index
    for cls in classes:
        if cls in pivot.columns:
            plt.plot(bins, pivot[cls], label=cls, color=colors[TE_CLASSES.index(cls)], linewidth=2)
    plt.title("Repetitive Proportion by Genetic Divergence/Class", fontsize=14)
    plt.xlabel("Genetic Divergence (%)")
    plt.ylabel("Proportion of Genome")
    # Adjusted x-axis limits for reduced offset
    plt.xlim(0, max_divergence + spacing / 2)
    plt.legend()
    plt.savefig(f"{output_prefix}_line.png", bbox_inches="tight")
    plt.close()

def create_pie_chart(df, classes, colors, genome_size, output_prefix):
    """Generate a pie chart."""
    total_proportions = df.groupby("Class")["Length"].sum()
    proportions = [total_proportions.get(cls, 0) / genome_size for cls in classes]
    unmasked = 1 - sum(proportions)
    proportions.append(unmasked)
    labels = classes + ["Unmasked"]
    colors.append("#000000")
    
    # Output the data to a .tsv file
    pie_data = pd.DataFrame({
        'Class': labels,
        'Proportion': proportions
    })
    output_file = f"{output_prefix}_pie_data.tsv"
    pie_data.to_csv(output_file, sep='\t', index=False)
    print(f"Data has been written to {output_file}")

    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(
        proportions, 
        labels=None, 
        colors=colors, 
        autopct='%1.1f%%', 
        startangle=140, 
        textprops={'color': "w"},
        wedgeprops={'edgecolor': 'white'}
    )
    
    # Add a legend
    ax.legend(
        wedges,
        labels,
        title="TE Classes",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1)
    )

    plt.title("Proportion of Genome Occupied by TE Classes", fontsize=14)
    plt.savefig(f"{output_prefix}_pie.png", bbox_inches="tight")
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Generate TE class plots from RepeatMasker output.")
    parser.add_argument("-r", "--repeatmasker", required=True, help="Compressed RepeatMasker .out.gz file")
    parser.add_argument("-s", "--summary", required=True, help="Summary .gz file to calculate genome size")
    parser.add_argument("--bin", type=int, default=1, help="Bin size for genetic divergence (default: 1)")
    parser.add_argument("--spacing", type=int, default=8, help="Spacing for x-axis (default: 8)")
    parser.add_argument("-c", "--classes", default="DIRS,LINE,SINE,LTR,RC,DNA", 
                        help="TE classes to include in plots, comma-separated (default: DIRS,LINE,SINE,LTR,RC,DNA)")
    parser.add_argument("-maxd", "--maxdivergence", type=int, default=40, 
                        help="Maximum genetic divergence for plots (default: 40)")
    parser.add_argument("-proc", "--processors", type=int, default=1, 
                        help="Number of processors for parallel processing (default: 1)")
    parser.add_argument("-op", "--output_prefix", required=True, help="Output prefix for plot filenames")
    args = parser.parse_args()

    classes = args.classes.split(",")
    genome_size = parse_summary_file(args.summary)
    df = parse_repeatmasker_out(args.repeatmasker, genome_size, args.bin, classes)
    
    create_stacked_bar_plot(df, classes, COLORS, args.maxdivergence, args.spacing, args.output_prefix)
    create_line_plot(df, classes, COLORS, args.maxdivergence, args.spacing, args.output_prefix)
    create_pie_chart(df, classes, COLORS[:len(classes)], genome_size, args.output_prefix)

if __name__ == "__main__":
    main()
