import argparse
import pickle

import matplotlib.pyplot as plt
import pandas as pd


def args():
    parser = argparse.ArgumentParser(description='Filters Blat output')
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument('-input_exonstat',
                               type=str,
                               help='path to exon stat pickle file input',
                               required=True)

    requiredNamed.add_argument(
        '-input_denovoexons',
        type=str,
        help='Path to tsv of denovo versus human database',
        required=True)

    requiredNamed.add_argument(
        '-input_denovo_nonexons',
        type=str,
        default=4,
        help='Path to tsv of denovo versus human alternative frame database',
        required=True)

    requiredNamed.add_argument('-output_pdf',
                               type=str,
                               help='Path to pdf plots output file',
                               required=True)

    return parser


cols = [
    "matches", "misMatches", "repMatches", "nCount", "qNumInsert",
    "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize",
    "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount",
    "blockSizes", "qStarts", "tStarts", "qseq", "tseq"
]

if __name__ == "__main__":
    arguments = args().parse_args()
    input_dn_exons = arguments.input_denovoexons
    input_dn_introns = arguments.input_denovo_nonexons
    input_exon_stat = arguments.input_exonstat
    output_pdf = arguments.output_pdf

    # loading exon stats
    exon_stat = pd.Series(pickle.load(open(input_exon_stat, "rb")))

    # loading dataframes
    df_dn_exons_td = pd.read_csv(input_dn_exons, sep="\t", header=0)

    df_dn_introns = pd.read_csv(input_dn_introns, sep="\t", header=0)

    df_dn_exons_td["n_aa"] = df_dn_exons_td.denovo_seq_nomod.map(
        lambda x: len(set(list(x))))
    df_dn_introns["n_aa"] = df_dn_introns.denovo_seq_nomod.map(
        lambda x: len(set(list(x))))

    df_dn_exons = df_dn_exons_td[~df_dn_exons_td.isDecoy]
    df_dn_exons_decoys = df_dn_exons_td[df_dn_exons_td.isDecoy]

    # plot 1 data prep
    data1 = [df_dn_exons.qName.nunique(), df_dn_introns.qName.nunique()]
    index = [
        f"Exonic spectra\n{df_dn_exons.qName.nunique()}",
        f"3FT spectra\n{df_dn_introns.qName.nunique()}"
    ]
    data1 = pd.Series(data1, index=index)

    # plot 2 data prep
    data2 = pd.DataFrame(
        {"Exonic spectra": df_dn_exons.predicted_score.tolist()})
    temp = df_dn_introns["predicted_score"].tolist()
    data2.loc[:len(temp) - 1, "3FT spectra"] = temp

    # plot 5
    data_scores = pd.DataFrame()
    data_scores["Target"] = df_dn_exons.predicted_score.tolist()
    temp = df_dn_exons_decoys.predicted_score.tolist()
    data_scores.loc[:len(temp) - 1, "Decoy"] = temp
    # plot 6
    n_aas = [1, 2, 3, 4]
    fdr = []
    index = []
    for n_aa in n_aas:
        ntargets = df_dn_exons[df_dn_exons.n_aa >= n_aa].shape[0]
        ndecoys = df_dn_exons_decoys[df_dn_exons_decoys.n_aa >= n_aa].shape[0]
        fdr.append(ndecoys / (ntargets + ndecoys))
        index.append(f"FDR {n_aa}+ unique amino acids")

    data4 = pd.Series(fdr, index=index)

    # plot 4 data prep
    data6 = pd.concat([
        df_dn_exons.denovo_seq_nomod.map(lambda x: len(set(x))).value_counts(),
        df_dn_exons_decoys.denovo_seq_nomod.map(
            lambda x: len(set(x))).value_counts()
    ],
                      axis=1)
    data6 = data6.fillna(0).astype(int)
    data6 = data6.apply(lambda x: x / x.sum())
    data6.columns = ["Target", "Decoy"]

    fig = plt.figure(figsize=(9, 6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3)

    # Denovo spectra piechart
    ax = fig.add_subplot(gs[0, 0])

    radius = 1.3
    width = 0.3
    # First Ring (outside)
    ax.axis('equal')
    ax.pie(data1.values,
           radius=radius,
           autopct='%1.1f%%',
           labels=data1.index.tolist(),
           colors=["blue", "black"],
           wedgeprops=dict(width=width, edgecolor='w'))

    # denovo spectra scores distribution #
    ax = fig.add_subplot(gs[0, 1])

    data2.plot(kind="hist",
               bins=20,
               density=1,
               ax=ax,
               color=["blue", "black"],
               legend=True,
               alpha=0.75,
               edgecolor='white')
    ax.set_title("Denovo spectra scores distribution")
    ax.set_xlabel("Denovo score")
    ax.set_ylabel("Fraction (probability density)")

    # 3FT spectra mismatch stats
    ax = fig.add_subplot(gs[0, 2])
    ax.set_title("3FT spectra mismatch stats", fontsize=8)
    exon_stat_plot_df = exon_stat.value_counts().sort_index()
    exon_stat_plot_df[exon_stat_plot_df.index >= 0].plot(kind="bar", ax=ax)
    ax.set_xlabel("3FT matches - Exon matches")
    ax.set_ylabel("Log(count)")
    ax.set_yscale("log")

    #  decoy vs targets scores plot
    ax = fig.add_subplot(gs[1, 0])
    data_scores.plot(kind="hist",
                     bins=40,
                     ax=ax,
                     color=["blue", "red"],
                     legend=True,
                     alpha=0.75,
                     edgecolor='white')
    ax.set_title("target decoy (1+ unique amino acid)\nscores distribution")
    ax.set_xlabel("Denovo predicted score")
    ax.set_ylabel("Frequency")

    # number of unique amino acids per spectrum plot
    ax = fig.add_subplot(gs[1, 1])
    data6.plot(kind="bar",
               ax=ax,
               color=["blue", "red"],
               legend=True,
               alpha=0.75)
    ax.set_title("Number of unique aa\nper spectrum distribution")
    ax.set_xlabel("Unique amino acids")
    ax.set_ylabel("Fraction of spectra")

    # fdr plot target decoy
    ax = fig.add_subplot(gs[1, 2])
    data4.plot(kind='bar', alpha=0.75, legend=False, ax=ax)
    ax.set_title("Target Decoy FDR")
    ax.set_xlabel("Denovo predicted score")
    ax.set_ylabel("FDR %\n#Decoy / (#Target + #Decoy)")
    ax.set_xticklabels([str(x) + '+' for x in n_aas])
    ax.set_xlabel('Number of unique amino acids per spectrum')

    for i, v in enumerate(data4.tolist()):
        ax.text(i,
                v + data4.max() * 0.01,
                f"{v*100:.2f}%",
                color='black',
                ha="center")

    plt.savefig(output_pdf, transparent=True)
