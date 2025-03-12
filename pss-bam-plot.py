#!/usr/bin/python

import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

params = { "figure.dpi": 500,
           "axes.titlesize": 25,
           "xtick.labelsize": 15,
           "ytick.labelsize": 15 }
plt.rcParams.update(params)     

nt_pairs = [ "AA", "AC", "AG", "AT",
             "CA", "CC", "CG", "CT",
             "GA", "GC", "GG", "GT",
             "TA", "TC", "TG", "TT" ]
sub_pairs = [p for p in nt_pairs if p not in ["AA", "CC", "GG", "TT"]]
color_dict = {"A": "#7bc043", "C": "#44a0f3", "G": "#ffd700", "T": "#db3401",
              "TC": "#8b0000", "AG": "#2a670f"} 

def load_args():
    '''Retrieves command line arguments'''
    desc = "pss-bam-plot.py: Create DNA damage plot (nucleotide composition & substitution) from pss-bam's output"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("-p", "--prefix", action="store", dest="input", required=True, metavar="STR", help=".pss.counts.txt and .pss.rates.txt file prefix (same as given to -o argument in pss-bam), output plot will also have this prefix")
    parser.add_argument("-r", "--region-length", action="store", dest="length", default=15, type=int, metavar="INT", help="length in basepairs into the interior of alignments to report (same as given to -r argument in pss-bam; default=15)")
    parser.add_argument("-m", "--max-rate", action="store", dest="rate", default=0.1, type=float, metavar="FLOAT", help="maximum substitution rate to set for y axis (default=0.1)")
    args = parser.parse_args()
    return args


def get_counts(in_prefix, region_len):
    fn = in_prefix + ".pss.counts.txt"
    f = open(fn, "r")
    n_skip = 0
    line = f.readline()
    while line and not line.startswith("### Reverse"):
        n_skip += 1
        line = f.readline()
    n_skip += 1
    
    fp_df = pd.read_table(fn, sep="\s+", comment="#", names=nt_pairs, nrows=region_len+2)
    tp_df = pd.read_table(fn, sep="\s+", skiprows=n_skip, names=nt_pairs, nrows=region_len+2)
    tp_df.index = np.arange(region_len-1, -3, -1)
    for df in [fp_df, tp_df]:
        df["A"] = df["AA"] + df["AC"] + df["AG"] + df["AT"]
        df["C"] = df["CA"] + df["CC"] + df["CG"] + df["CT"]
        df["G"] = df["GA"] + df["GC"] + df["GG"] + df["GT"]
        df["T"] = df["TA"] + df["TC"] + df["TG"] + df["TT"]
    f.close()
    return fp_df, tp_df


def get_rates(in_prefix, region_len):
    fn = in_prefix + ".pss.rates.txt"
    f = open(fn, "r")
    n_skip = 0
    line = f.readline()
    while line and not line.startswith("### Reverse"):
        n_skip += 1
        line = f.readline()
    n_skip += 1
    
    fp_df = pd.read_table(fn, sep="\s+", comment="#", names=sub_pairs, nrows=region_len, dtype=float)
    tp_df = pd.read_table(fn, sep="\s+", skiprows=n_skip, names=sub_pairs, nrows=region_len, dtype=float)
    tp_df.index = np.arange(region_len-1, -1, -1)
    f.close()
    return fp_df, tp_df


def make_plot(fp_counts, tp_counts, fp_rates, tp_rates, in_prefix, region_len, max_rate):
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(12, 8))
    plt.subplots_adjust(wspace=0.15)
    
    for ax in [ax1, ax2]:
        ax.set_xlim(-3, region_len)
        ax.set_xticks(np.arange(-2, region_len, 1))
        ax.set_ylim(0, max_rate)
        ax.add_patch(mpatches.Rectangle((-3, 0), 2.5, max_rate, color="silver", zorder=0))
    ax1.set_xticklabels(np.arange(-2, region_len, 1), fontsize=13)
    ax1.set_ylabel("Substitution rate", labelpad=10, fontsize=20)
    ax1.set_title("5' end", pad=8)

    ax2.set_xticklabels(np.concatenate((np.array([2,1]), np.arange(0, region_len, 1))), fontsize=13)
    ax2.invert_xaxis()
    ax2.set_title("3' end", pad=8)
    
    for x in range(-2, region_len):
        fp_scale = max_rate/(fp_counts.at[x,"A"] + fp_counts.at[x,"C"] + fp_counts.at[x,"G"] + fp_counts.at[x,"T"])
        tp_scale = max_rate/(tp_counts.at[x,"A"] + tp_counts.at[x,"C"] + tp_counts.at[x,"G"] + tp_counts.at[x,"T"])
        fy = 0
        ty = 0
        for b in ["A", "G", "C", "T"]:
            ax1.bar(x, fp_counts.at[x,b]*fp_scale, bottom=fy, color=color_dict[b], edgecolor="black")
            fy += fp_counts.at[x,b]*fp_scale
            ax2.bar(x, tp_counts.at[x,b]*tp_scale, bottom=ty, color=color_dict[b], edgecolor="black")
            ty += tp_counts.at[x,b]*tp_scale

    for p in sub_pairs:
        if p == "TC":
            line1, = ax1.plot(fp_rates[p], color=color_dict[p], lw=3, label="C>T")
            ax2.plot(tp_rates[p], color=color_dict[p], lw=3)
        elif p == "AG":
            line2, = ax1.plot(fp_rates[p], color=color_dict[p], lw=3, label="G>A")
            ax2.plot(tp_rates[p], color=color_dict[p], lw=3)     
        else:
            line3, = ax1.plot(fp_rates[p], color="black", lw=0.75, label="Others")
            ax2.plot(tp_rates[p], color="black", lw=0.75)

    legend_handles = [line1, line2, line3]
    for b in ["A", "G", "C", "T"]:
        patch = mpatches.Patch(color=color_dict[b], label=b)
        legend_handles.append(patch)
    
    plt.legend(handles=legend_handles, bbox_to_anchor=(0.95,-0.05), ncol=8, frameon=False, fontsize=18)
    
    out_fn = in_prefix + ".pss.plot.svg"
    plt.savefig(fname=out_fn, format="svg", dpi=500)
    return None


def main():
    args = load_args()
    in_prefix = args.input
    region_len = args.length
    max_rate = args.rate
    fp_counts, tp_counts = get_counts(in_prefix, region_len)
    fp_rates, tp_rates = get_rates(in_prefix, region_len)
    make_plot(fp_counts, tp_counts, fp_rates, tp_rates, in_prefix, region_len, max_rate)


if __name__ == "__main__":
    main()
