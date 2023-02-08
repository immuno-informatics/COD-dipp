import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.patches import Rectangle

import argparse
import os


def args():
    parser = argparse.ArgumentParser(
        description=
        "Generates plot after TITER analysis(doi: 10.1093/bioinformatics/btx247)"
    )
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument('-input_titer_prediction',
                               type=str,
                               help='Path to tsv of TITER output',
                               required=True)

    requiredNamed.add_argument(
        '-input_frames',
        type=str,
        help='Path to tsv with positive strand features',
        required=True)

    requiredNamed.add_argument('-output',
                               type=str,
                               help='Path to output pdf file',
                               required=True)

    return parser


def write_panellabels(axes):
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ".lower()
    for i, ax in enumerate(axes):
        label = alphabet[i] + ")"
        ax.text(-0.1,
                1.15,
                label,
                transform=ax.transAxes,
                fontsize=8,
                fontweight='bold',
                va='top',
                ha='right')


CODONS = {"ATG": "M", "CTG": "L", "GTG": "V", "ACG": "T", "TTG": "L"}


def sunburst(nodes, total=np.pi * 2, offset=0, level=0, ax=None):
    ax = ax or plt.subplot(111, projection='polar')

    if level == 0 and len(nodes) == 1:
        label, value, subnodes = nodes[0]
        ax.bar([0], [0.5], [np.pi * 2])
        ax.text(0, 0, label, ha='center', va='center')
        sunburst(subnodes, total=value, level=level + 1, ax=ax)
    elif nodes:
        d = np.pi * 2 / total
        labels = []
        widths = []
        local_offset = offset
        for label, value, subnodes in nodes:
            labels.append(label)
            widths.append(value * d)
            sunburst(subnodes,
                     total=total,
                     offset=local_offset,
                     level=level + 1,
                     ax=ax)
            local_offset += value
        values = np.cumsum([offset * d] + widths[:-1])
        heights = [1] * len(nodes)
        bottoms = np.zeros(len(nodes)) + level - 0.5
        rects = ax.bar(values,
                       heights,
                       widths,
                       bottoms,
                       linewidth=1,
                       edgecolor='white',
                       align='edge')
        for rect, label in zip(rects, labels):
            x = rect.get_x() + rect.get_width() / 2
            y = rect.get_y() + rect.get_height() / 2
            rotation = (90 + (360 - np.degrees(x) % 180)) % 360
            ax.text(x, y, label, rotation=rotation, ha='center', va='center')

    if level == 0:
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location('N')
        ax.set_axis_off()
    return ax


if __name__ == "__main__":
    arguments = args().parse_args()
    file_raw = arguments.input_titer_prediction
    input_frames = arguments.input_frames
    output = arguments.output

    df_raw = pd.read_csv(file_raw, sep="\t")
    df_introns = pd.read_csv(input_frames, sep="\t")

    df_raw["chr_coords"] = df_raw["chr_coords"].map(str)
    df_introns["chr_coords"] = df_introns["chr_coords"].map(str)

    df = df_raw[df_raw["pass threshold"]].sort_values(
        'titer score',
        ascending=False).drop_duplicates("denovo_peptide_coords").reset_index(
            drop=True)

    if df.empty:
        print(
            "No intronic peptide passes the TIS (Translation Initiation site threshold). exiting ..."
        )
        os.system(f"touch {output}")
        exit(0)

    df_raw_unique = df_raw.sort_values("titer score",
                                       ascending=False).drop_duplicates([
                                           "chr_coords", "start_coords",
                                           "end_coords", "strand_coords"
                                       ])
    df_introns = df_introns.merge(df_raw_unique, how="left")
    df_introns["with start codon"] = ~df_introns["titer seq"].isnull()
    df_introns.columns = df_introns.columns.str.replace(
        "pass threshold", "with tis")
    df_introns["with tis"] = df_introns["with tis"].fillna(False)
    df_introns["start codon"] = df_introns["titer seq"].map(
        lambda x: x[100:103] + f" ({CODONS[x[100:103]]})"
        if isinstance(x, str) else "NULL")

    inframe_count = df_introns["inframe intron"].value_counts(
        normalize=True)[True] * 100
    with_start_codon = df_introns["with start codon"].value_counts(
        normalize=True)[True] * 100
    start_codons = df_introns["start codon"].value_counts(normalize=True) * 100
    with_tis = df_introns.groupby("start codon").agg({
        "with tis": "value_counts"
    }).unstack("start codon").loc[True, :].reset_index(level=0, drop=True)
    with_tis = (with_tis / df_introns.shape[0]) * 100

    # yapf: disable
    # ('Inframe', inframe_count, []),
    data = [
        ('Intronic\npeptides', 100, [
            (f"With upstream\nstart codon\n{with_start_codon:.2f}%", with_start_codon,
                [
                    (f"ATG (M)\n{start_codons['ATG (M)']:.2f}%", start_codons["ATG (M)"],
                    [
                        (f"TIS\n{with_tis['ATG (M)']:.2f}%", with_tis["ATG (M)"], [])
                    ]),
                    (f"CTG (L)\n{start_codons['CTG (L)']:.2f}%", start_codons["CTG (L)"],
                     [
                        (f"TIS\n{with_tis['CTG (L)']:.2f}%", with_tis["CTG (L)"], [])
                    ]),
                    (f"ACG (T)\n{start_codons['ACG (T)']:.2f}%", start_codons["ACG (T)"],
                     [
                        (f"TIS\n{with_tis['ACG (T)']:.2f}%", with_tis["ACG (T)"], [])
                    ]),
                    (f"GTG (V)\n{start_codons['GTG (V)']:.2f}%", start_codons["GTG (V)"],
                     [
                        (f"TIS\n{with_tis['GTG (V)']:.2f}%", with_tis["GTG (V)"], [])
                    ]),
                    (f"TTG (L)\n{start_codons['TTG (L)']:.2f}%", start_codons["TTG (L)"],
                     [
                        (f"TIS\n{with_tis['TTG (L)']:.2f}%", with_tis["TTG (L)"], [])
                    ]),
                ]),
        ])
    ]

    # yapf: enable
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), constrained_layout=True)
    ax = sunburst(data, ax)
    plt.savefig(output, transparent=True)

# fig, axes = plt.subplots(2, 3, figsize=(7, 5), constrained_layout=True)
# fig.suptitle("TITER analysis of HLA-I denovo intronic peptides")

# axes = axes.flatten()
# write_panellabels(axes)

# ax = axes[0]
# df_introns = df_introns.sort_values("inframe intron", ascending=False)
# inframe = df_introns.drop_duplicates(
#     "denovo_peptide_coords")["inframe intron"].value_counts()
# inframe_ratio = df_introns.drop_duplicates(
#     "denovo_peptide_coords")["inframe intron"].value_counts(normalize=True)
# inframe_ratio.plot(kind="bar", ax=ax)
# for x, value in enumerate(inframe_ratio):
#     ax.text(x, value + 0.025, inframe.iloc[x], fontsize=8, ha='center')

# ax.set_xlabel("Inframe with upstream exon")
# ax.set_ylabel("Fraction of unique intronic peptides")
# ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
# lg_title = "Inframe with\nupstream exon"

# ax = axes[1]
# df_introns = df_introns.sort_values("with start codon", ascending=False)
# with_tis = df_introns.drop_duplicates("denovo_peptide_coords")[[
#     "with start codon", "inframe intron"
# ]].value_counts().unstack()
# with_tis_ratio = df_introns.drop_duplicates("denovo_peptide_coords")[[
#     "with start codon", "inframe intron"
# ]].value_counts(normalize=True).unstack()
# with_tis_ratio.plot(kind="bar", ax=ax, stacked=True, legend=True)
# for x, value in enumerate(with_tis_ratio.sum(axis=1)):
#     ax.text(x,
#             value + 0.025,
#             with_tis.iloc[x].sum(),
#             fontsize=8,
#             ha='center')

# ax.set_xlabel("Intronic peptide with upstream start codon")
# ax.set_ylabel("Fraction of unique intronic peptides")
# ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
# ax.legend(ax.get_legend_handles_labels()[0],
#           ax.get_legend_handles_labels()[1],
#           title=lg_title,
#           title_fontsize=8,
#           fontsize=8)

# ax = axes[2]
# # histogram score
# df_scores = pd.concat([
#     df_raw[df_raw["pass threshold"]]["titer score"],
#     df_raw[df_raw["pass threshold"] == False]["titer score"]
# ],
#                       axis=1,
#                       ignore_index=True)
# df_scores.columns = ["Pass threshold = True", "Pass threshold = False"]
# bins = np.arange(df_scores.min().min(), df_scores.max().max(), 0.2)
# colors = ["b", "r"]
# ax.hist(df_scores["Pass threshold = True"],
#         bins=bins,
#         color=colors[0],
#         edgecolor="white")
# ax.hist(df_scores["Pass threshold = False"],
#         bins=bins,
#         color=colors[1],
#         edgecolor="white")
# ax.set_xlabel("TIS prediction score (TITER tool)")
# ax.set_ylabel("Frequency")

# #create legend
# handles = [Rectangle((0, 0), 1, 1, color=c) for c in colors]
# ax.legend(handles, df_scores.columns)

# ax = axes[3]

# df_introns = df_introns.sort_values("with tis", ascending=False)

# tis_df = df_introns[df_introns["with start codon"]].drop_duplicates(
#     "denovo_peptide_coords")[["with tis",
#                               "inframe intron"]].value_counts().unstack()
# tis_df_ratio = df_introns[df_introns["with start codon"]].drop_duplicates(
#     "denovo_peptide_coords")[["with tis", "inframe intron"
#                               ]].value_counts(normalize=True).unstack()
# tis_df_ratio.plot(kind="bar", ax=ax, stacked=True, legend=True)
# for x, value in enumerate(tis_df_ratio.sum(axis=1)):
#     ax.text(x,
#             value + 0.025,
#             tis_df.iloc[x].sum(),
#             fontsize=8,
#             ha='center')

# ax.set_xlabel("Intronic peptide with upstream TIS")
# ax.set_ylabel("Frequency")
# ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
# ax.legend(ax.get_legend_handles_labels()[0],
#           ax.get_legend_handles_labels()[1],
#           title=lg_title,
#           title_fontsize=8,
#           fontsize=8)

# ax = axes[4]
# df["titer seq"].map(
#     lambda x: x[100:103] + f" ({CODONS[x[100:103]]})").value_counts().plot(
#         kind="bar", ax=ax, color="b")
# ax.set_xlabel("TIS start codon")
# ax.set_ylabel("Frequency")

# plt.savefig(output, transparent=True)
