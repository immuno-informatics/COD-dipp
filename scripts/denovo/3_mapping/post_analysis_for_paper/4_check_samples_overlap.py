"""
    Checks if Denovo intronic peptides are shared between samples
    since authors DOI: 10.1126/scitranslmed.aau5516 claimed that
    TSAs (Tumor Specific Antigens) comming from non coding regiosn are shared
    between patients
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as ticker
from math import log
from scipy.cluster import hierarchy
from meta_variables import COLOR_CYCLER, DISEASE_ORDER

NETMHCPAN = "/net/archive/groups/plgg_iccvs/tools/netmhccons/netMHCpan-4.0/netMHCpan"
NETMHCPAN_PARSER = "$HOME/git_repositories/covid19/src/NetMHCPan/netmhcpan_parser.py"
POPULATION_COVERAGE_SCRIPT = "/net/archive/groups/plgg_iccvs/tools/population_coverage/calculate_population_coverage.py"
IMMUNOPRED = "/net/archive/groups/plgg_iccvs/tools/immunogenicity/predict_immunogenicity.py"


def filter_peptides(x):
    """
        keeps peptides detected at least 2 times in a sample
        returns the peptide count per disease type
    """
    npsample = x["Sample Name"].value_counts()
    index = npsample[npsample >= 5].index
    subx = x[x["Sample Name"].isin(index)]
    return subx.drop_duplicates("Sample Name")["Disease"].value_counts()


def filter_peptides_bysample(x):
    """
        keeps peptides detected at least 2 times in a sample
        returns the peptide count per disease type
    """
    npsample = x["Sample Name"].value_counts()
    index = npsample[npsample >= 5].index
    subx = x[x["Sample Name"].isin(index)]
    return subx["Sample Name"].value_counts()


def heatmap_of_recurrent_pancancer_cryptic_peptides_depricated(
        best_peptides, count_dict, output_pdf):
    plot_df = best_peptides * 100
    plot_df = plot_df[~plot_df.index.isin(["Na"])]

    # min_sample_count = 5
    # count_series = pd.Series(count_dict)
    # plot_df = plot_df[plot_df.index.isin(count_series[count_series >= min_sample_count].index)]

    # index = [
    #     "Disease Free", "Colon Cancer",
    #     "Metastatic Malignant Melanoma", "Melanoma", "Hepatocellular Carcinoma",
    # ]
    #
    # plot_df = best_peptides.loc[best_peptides.index.isin(index), :] * 100

    plot_df = plot_df.loc[:, plot_df.loc["Disease Free", :] == 0]
    g = sns.clustermap(plot_df, figsize=(11, 6), cmap="YlGnBu")
    g.fig.axes[0].clear()
    ydata = pd.Series(count_dict).loc[[
        x.get_text() for x in g.ax_heatmap.get_yticklabels()
    ]].tolist()
    y_pos = range(len(ydata))
    ax = g.fig.axes[0]
    ax.barh(y_pos, ydata, height=-0.5)

    ax.invert_yaxis()  # labels read top-to-bottom
    ax.invert_xaxis()  # labels read top-to-bottom

    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks([])
    ax2.set_yticklabels([])

    ax.axes.get_yaxis().set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax2.set_xlabel("Number\nof samples")
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), fontsize=6)
    plt.savefig(output_pdf, bbox_tight=True)
    return plot_df, g


def color_generation(df_covariates, alpha=0.75):
    color_dict = {}
    values = pd.Series(df_covariates.values.flatten()).drop_duplicates().values
    for index, item in enumerate(values):
        x = COLOR_CYCLER[index]
        color_dict[item] = (x[0], x[1], x[2], alpha)
    return df_covariates.applymap(lambda x: color_dict.get(x, np.nan))


def cluster_samples(df, groupby, clusterby):
    df_list = []
    for group, sdf in df.groupby(groupby):
        try:
            Z = hierarchy.linkage(sdf[clusterby].values)
            df_list.append(sdf.iloc[list(hierarchy.leaves_list(Z))])
        except ValueError:
            df_list.append(sdf)
    return pd.concat(df_list)


def add_spaces(df, space, groupby, order=None):
    spacer = pd.DataFrame(index=["spacer"] * space,
                          columns=df.columns).fillna(0)
    df_dict = {}
    for group, sdf in df.groupby(groupby):
        spacer[groupby] = group
        df_dict[group] = pd.concat([sdf, spacer])
    df_list = []
    if order:
        for group in order:
            try:
                df_list.append(df_dict[group])
            except KeyError:
                continue
    else:
        df_list = list(df_dict.values())
    return pd.concat(df_list)


def pathogenic_cryptic_peptides(df, df_ann, min_sample_recurrence):
    ann_cols = ["Sample Name", "Sample Type", "Disease", "Treatment"]
    df_pat = df.groupby("peptide").apply(lambda x: x.drop_duplicates(
        "Sample Name").set_index("Sample Name")["count"]).unstack(
            level=1).T.fillna(0)
    df_pat = df_pat.reset_index().merge(df_ann[ann_cols], how="left")
    peptides = df_pat.columns.tolist()[1:-3]
    temp = df_pat[df_pat["Disease"] == "Disease Free"][peptides].sum()
    df_peptides = df_pat.drop(temp[temp != 0].index.tolist(), axis=1)
    df_peptides = df_peptides[df_peptides["Disease"] != "Disease Free"]
    peptides = df_peptides.columns.tolist()[1:-3]
    return df_peptides[peptides].groupby(
        df_peptides.Disease).apply(lambda x: x.sum()).T


def heatmap_of_recurrent_pancancer_cryptic_peptides(df, df_ann,
                                                    min_sample_recurrence,
                                                    count_dict, output_pdf):
    pdf = PdfPages(output_pdf)
    ann_cols = [
        "Sample Name", "Sample Type", "Disease", "Treatment", "Pride ID"
    ]
    plot_df = df.groupby("peptide").apply(lambda x: x.drop_duplicates(
        "Sample Name").set_index("Sample Name")["count"]).unstack(
            level=1).T.fillna(0)
    plot_df.index.name = "Sample Name"
    heatmap_cols = plot_df.columns.tolist()
    plot_df = plot_df.reset_index().merge(df_ann[ann_cols], how="left")
    df_col = color_generation(plot_df[ann_cols[1:]])
    df_col.loc["spacer", :] = np.nan
    df_col.loc["spacer", :] = df_col.loc["spacer", :].map(
        lambda x: (1.0, 1.0, 1.0) if x is np.nan else x)
    plot_df = cluster_samples(plot_df, "Disease", heatmap_cols)
    plot_df = add_spaces(plot_df, 5, "Disease", order=DISEASE_ORDER)
    index_heatmap = plot_df.index != "spacer"
    plot_df.loc[index_heatmap, heatmap_cols] = plot_df.loc[
        plot_df.index != "spacer", heatmap_cols].applymap(lambda x: log(x + 1))
    g = sns.clustermap(plot_df[heatmap_cols],
                       figsize=(9, 12),
                       cmap="binary",
                       row_colors=df_col,
                       row_cluster=False)
    g.ax_heatmap.set_yticklabels([])
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.set_xlabel("Recurrent cryptic peptides")
    disease_count_ordered = plot_df["Disease"].value_counts().loc[
        plot_df.Disease.unique()]
    ticks_major = disease_count_ordered.cumsum()
    ticks_minor = ticks_major - (disease_count_ordered / 2)
    g.ax_heatmap.set_yticks(ticks_major.tolist())
    g.ax_heatmap.yaxis.set_major_formatter(ticker.NullFormatter())
    g.ax_heatmap.yaxis.set_minor_locator(
        ticker.FixedLocator(ticks_minor.tolist()))
    g.ax_heatmap.tick_params(axis=u'both', which=u'both', length=0)
    g.ax_heatmap.yaxis.set_minor_formatter(
        ticker.FixedFormatter(ticks_minor.index.tolist()))
    # yapf: disable
    [g.ax_heatmap.hlines(y - 2.5, 0, plot_df.shape[1], colors="black", linewidths=0.5) for y in ticks_major.values]
    # yapf: enable
    g.ax_col_dendrogram.set_visible(False)
    pdf.savefig(bbox_tight=True)
    fig, ax = plt.subplots(1, 1, figsize=(6, 9), constrained_layout=True)
    # add legends
    lgs = []
    for i, col in enumerate(ann_cols[1:]):
        titles = plot_df.loc[index_heatmap, col].drop_duplicates()
        values = df_col.loc[titles.index, col].drop_duplicates()
        handles = [Patch(facecolor=value) for value in values.tolist()]
        lgs.append(
            ax.legend(handles,
                      titles.tolist(),
                      loc=i + 1,
                      title=col,
                      fontsize=8))
    [ax.add_artist(lg) for lg in lgs[:-1]]
    pdf.savefig(bbox_tight=True)
    pdf.close()
    return plot_df, g


def get_hla_binding(peptide_list, supertypes=None):
    pep_file = "/tmp/pepfile.txt"
    netmhcpan_output = "/tmp/pepfile_netmhcpan"
    netmhcpan_parsed = "/tmp/pepfile_netmhcpan_parsed"

    if not supertypes:
        supertypes = "HLA-A02:01,HLA-A01:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01,HLA-B15:01"

    with open(pep_file, "w") as fh:
        fh.write("\n".join(peptide_list))

    os.system(
        f"{NETMHCPAN} -f {pep_file} -a {supertypes} -inptype 1 -l 8,9,10,11,12 > {netmhcpan_output}"
    )
    os.system(
        f"python {NETMHCPAN_PARSER} {netmhcpan_output} {netmhcpan_parsed}")

    bp_binding = pd.read_csv(netmhcpan_parsed, sep="\t")
    os.system(f"rm -rf {pep_file} {netmhcpan_output} {netmhcpan_parsed}")
    return bp_binding


def population_coverage_plot(bp_binding, output_dir, output_log):
    pep_file = "/tmp/pepfile.txt"
    allele_coverage_input = bp_binding[~bp_binding.BindLevel.isnull()]
    allele_coverage_input = allele_coverage_input.groupby("Peptide").apply(
        lambda x: ",".join(x.HLA.tolist()))
    allele_coverage_input.to_csv(pep_file, index=True, header=False, sep="\t")

    os.system(
        f"python {POPULATION_COVERAGE_SCRIPT} -f {pep_file} -p World -c I --plot {output_dir} > {output_log} 2>&1"
    )


def get_immunogenicity(peptides):
    """
        calculates_immunogenicity scores ref. Calis et. al. 2013
        input: (list) of peptides equences
        returns: (pd.DataFrame)
    """
    in_tmp = "/tmp/immunopred.txt"
    out_tmp = "/tmp/immunopred_res.txt"
    with open(in_tmp, "w") as out:
        out.write("\n".join(peptides))
    os.system(f"python {IMMUNOPRED} {in_tmp} > {out_tmp}")
    df = pd.read_csv(out_tmp, skiprows=3, header=0, sep=",")
    os.system(f"rm {in_tmp} {out_tmp}")
    return df


def plot_immunogenicity_scores(pep_list1, pep_list2, output_pdf):
    df_im1 = get_immunogenicity(pep_list1)
    df_im2 = get_immunogenicity(pep_list2)

    df_im1["tag"] = "all cryptic peptides"
    df_im2["tag"] = "Shared cryptic peptides"
    df_im = pd.concat([df_im1, df_im2], axis=0)

    fig, ax = plt.subplots(1, figsize=(4, 4), constrained_layout=True)
    for group, sdf in df_im.groupby("tag"):
        ax.hist(sdf.score,
                bins=10,
                label=group,
                density=True,
                edgecolor="white",
                alpha=0.75)

    ax.legend()
    plt.savefig(output_pdf, transparent=True)
    return df_im


def main():
    # yapf: disable
    root = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all"
    file_peptides = os.path.join(root, "peptides_probamconvert/gene_inference/bam_files_2/coverage/merged.sorted_features.tsv")
    # file_ann = os.path.join(root, "tables/sample_centered_table.tsv")
    file_ann = os.path.join(root, "tables/HLA-I_tables/samples_info.tsv")

    output_root = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new"
    output_pat_peptides = os.path.join(output_root, "pathogenic_peptides.tsv")
    output_bindings = os.path.join(output_root, "cryptic_peptides_bindings.tsv")
    output_bindings_best = os.path.join(output_root, "bestcryptic_peptides_bindings.tsv")
    output_pdf = os.path.join(output_root, "bestcryptic_cancer_clustermap.pdf")
    population_coverage_output = os.path.join(output_root, "cryp_peptide_worldcoverage.txt")
    # yapf: enable

    df_all = pd.read_csv(file_peptides, sep="\t")
    df = df_all[df_all["tag"] == "denovo cryptic"].copy(deep=True)
    df.columns = df.columns.str.replace("sample_name", "Sample Name")
    df.drop("study_id", axis=1, inplace=True)

    df_ann = pd.read_csv(file_ann, sep="\t")
    df_ann = df_ann[df_ann["HLA Class"] == 1]
    df_ann.drop("Disease", axis=1, inplace=True)
    df_ann.columns = df_ann.columns.str.replace("^Disease name$", "Disease")

    df = df.merge(df_ann)

    # adding the number of samples per shared pepties
    min_sample_recurrence = 10
    min_psm_per_sample = 2
    shared_peptides = df.groupby("peptide").apply(lambda x: sum(
        x.drop_duplicates("Sample Name")["count"] >= min_psm_per_sample) >=
                                                  min_sample_recurrence)
    shared_peptides = shared_peptides[shared_peptides].index

    df = df[df.peptide.isin(shared_peptides)]

    # disease sample counts
    count_dict = df_ann["Disease"].value_counts().to_dict()

    # identifying cryptic recurrent peptides
    # best_peptides = df[df["Number of samples"] >= min_sample_recurrence]
    # best_peptides = best_peptides.groupby("peptide").apply(filter_peptides)
    # best_peptides = best_peptides / best_peptides.index.map(
    #     lambda x: count_dict[x[1]])
    # best_peptides = best_peptides.unstack(level=1).T.fillna(0)

    # plotting heatmap of pan cancer view
    # df_pan, heatmap_pan = heatmap_of_recurrent_pancancer_cryptic_peptides_depricated(
    #     best_peptides, count_dict, output_pdf)
    df_pan, heatmap_pan = heatmap_of_recurrent_pancancer_cryptic_peptides(
        df, df_ann, min_sample_recurrence, count_dict, output_pdf)

    df_pathogenic = pathogenic_cryptic_peptides(df, df_ann,
                                                min_sample_recurrence)
    df_pathogenic.reset_index().to_csv(output_pat_peptides,
                                       sep="\t",
                                       header=True,
                                       index=False)
    # get the recurrent cryptic peptides bindings affinities to HLA supertype
    bp_binding = get_hla_binding(list(df.peptide.unique()))
    bp_binding_best = bp_binding[~bp_binding["BindLevel"].isnull(
    )].sort_values("%Rank").drop_duplicates('Peptide')

    # write binding affinities to tsv files
    bp_binding.to_csv(output_bindings, sep="\t", header=True, index=False)
    bp_binding_best.to_csv(output_bindings_best,
                           sep="\t",
                           header=True,
                           index=False)

    # plot the population coverage graph of the recurrent cryptic peptides
    population_coverage_plot(bp_binding, output_root,
                             population_coverage_output)

    # plot immunogenicity score cryptic vs non cryptic peptides
    output_pdf = os.path.join(output_root, "bestcryptic_immuno_scores.pdf")
    all_cryptic_peptides = list(
        df_all[(df_all["tag"] == "denovo cryptic")
               & (df_all.peptide.str.len() == 9)].peptide.unique())
    shared_crpytic_peptides = list(
        df[df.peptide.str.len() == 9].peptide.unique())
    df_immuno = plot_immunogenicity_scores(all_cryptic_peptides,
                                           shared_crpytic_peptides, output_pdf)
    return df_pathogenic, df_pan, heatmap_pan, bp_binding_best, df_immuno


if __name__ == "__main__":
    df_pat, df_pan, heatmap_pan, bp_binding_best, df_immuno = main()
    # TODO:
    #     - descriptive analysis about the number of SB and WB (SB: 181, WB: 35, shared peptides: 239)
    #     - genes envolvement in the shared peptides
