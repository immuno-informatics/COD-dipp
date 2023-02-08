import argparse

import matplotlib.pyplot as plt
import pandas as pd


def args():
    parser = argparse.ArgumentParser(
        description="Plots denovo mapped peptides stats")
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-input_coordinates_annotation',
        type=str,
        nargs=2,
        help='Path to annotated denovo peptides TSV files (order: 3m, 4m)',
        required=True)

    requiredNamed.add_argument('-output_pdf',
                               type=str,
                               help='path to output pdf file',
                               required=True)

    return parser


if __name__ == "__main__":
    arguments = args().parse_args()
    coor_ann_3m = arguments.input_coordinates_annotation[0]
    coor_ann_4m = arguments.input_coordinates_annotation[1]
    output = arguments.output_pdf

    df_ann_3m = pd.read_csv(coor_ann_3m, sep="\t")
    df_ann_4m = pd.read_csv(coor_ann_4m, sep="\t")

    data_3m = df_ann_3m.groupby("id_coords").agg(
        {"feature_ann":
         lambda x: list(set(x))})["feature_ann"].explode().value_counts()
    data_4m = df_ann_4m.groupby("id_coords").agg(
        {"feature_ann":
         lambda x: list(set(x))})["feature_ann"].explode().value_counts()

    data_3m.index = data_3m.index.str.replace(".", "Junction spanning")
    data_3m.index = data_3m.index.str.replace("^Exon$", "Exon\n(out of frame)")
    data_4m.index = data_4m.index.str.replace(".", "Junction spanning")
    data_4m.index = data_4m.index.str.replace("^Exon$", "Exon\n(out of frame)")

    fig, axes = plt.subplots(2, 1, figsize=(3, 6), constrained_layout=True)
    fig.suptitle(
        "Denovo peptides 3FT features barplot with at least n\namino acid mismatches from a protein sequence",
        fontsize=8)

    ax = axes[0]
    data_3m.plot(kind="bar", ax=ax)
    ax.set_ylabel("Frequency")
    ax.set_title("3 amino acid mismatch")
    for i, v in enumerate(data_3m.tolist()):
        ax.text(i,
                v + data_3m.max() * 0.01,
                str(v),
                color='black',
                ha="center")

    ax = axes[1]
    data_4m.plot(kind="bar", ax=ax)
    ax.set_ylabel("Frequency")
    ax.set_title("4 amino acid mismatch")
    for i, v in enumerate(data_4m.tolist()):
        ax.text(i,
                v + data_4m.max() * 0.01,
                str(v),
                color='black',
                ha="center")

    plt.savefig(output, transparent=True)
