import os

import numpy as np
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


def gen_figure(nrows, ncols):
    fig, axes = plt.subplots(nrows, ncols, figsize=(
        7.9, 11), constrained_layout=True)
    axes = axes.flatten()
    return fig, axes


def read_fasta_header(fasta_file):
    dict_seq = {}
    for bioseq in SeqIO.parse(fasta_file, format="fasta"):
        if bioseq.name.startswith("ENSP"):
            header = bioseq.name.split("|")
            try:
                dict_seq[header[1].rsplit(".", 1)[0]] = bioseq.description
            except IndexError:
                continue
    return dict_seq


aa_conv = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Leu/Ile": "L"
}

# link: 10.6084/m9.figshare.16538097
file_anno = "opensearch_PTMiner_anno_HLA-I.tsv.gz"
# link: 10.6084/m9.figshare.16538097
file_fil = "opensearch_PTMiner_filtered_FLR_HLA-I.tsv.gz"
# link: 10.6084/m9.figshare.16538097
file_loc = "opensearch_PTMiner_loc_HLA-I.tsv.gz"
# link: 10.6084/m9.figshare.16538097 in pipeline_annotation_files.zip
database_file = "2019-04-30-td-Homo_sapiens_GRCh38_biomart.fasta"
file_samples = "sample_centered_table.tsv"  # in article supp data

df_fil = pd.read_csv(file_fil, sep="\t")
df_loc = pd.read_csv(file_loc, sep="\t")
df_ann = pd.read_csv(file_anno, sep="\t")
db_header = read_fasta_header(database_file)
df_samples = pd.read_csv(file_samples, sep="\t")
# file_pub_variants = "../publication_variants.tsv"
# df_pub_var = pd.read_csv(file_pub_variants, sep="\t")
# df_pub_var.columns = df_pub_var.columns.str.replace("variation", "Variation")
# df_pub_var.columns = df_pub_var.columns.str.replace("transcript",
#                                                     "Transcript Access")

HLAI_samples = df_samples[df_samples["HLA Class"] == 1]["Sample Name"].tolist()

df_fil = df_fil[df_fil["sample"].isin(HLAI_samples)].reset_index(
    drop=True).copy(deep=True)
df_loc = df_loc[df_loc["sample"].isin(HLAI_samples)].reset_index(
    drop=True).copy(deep=True)
df_ann = df_ann[df_ann["sample"].isin(HLAI_samples)].reset_index(
    drop=True).copy(deep=True)


df_ann = df_ann[
    df_ann["Annotated Mod Classification"] != "Artefact"].reset_index(
        drop=True).copy(deep=True)
df_ann_saav = df_ann[df_ann["Annotated Mod"].str.contains(
    "->", na=False)].reset_index(drop=True).copy(deep=True)

df_ann_saav["Ref AA"] = df_ann_saav["Annotated Mod"].map(
    lambda x: x.split("->")[0])
df_ann_saav["Alt AA"] = df_ann_saav["Annotated Mod"].map(
    lambda x: x.split("->")[1].split(" ")[0])
df_ann_saav["Variation"] = df_ann_saav["Ref AA"] + df_ann_saav[
    "Protein Position"].map(str) + df_ann_saav["Alt AA"]
df_ann_saav["Transcript Access"] = df_ann_saav["Protein Access"].map(
    lambda x: x.split("|")[1] if x.startswith("ENSP") else np.nan)


fig, ax = plt.subplots(1, 1, figsize=(1.5, 1.5), constrained_layout=True)
df_ann.drop_duplicates(['Sequence', 'sample', 'pride id', 'Annotated Mod'])["Annotation Type"]\
    .value_counts().plot(kind="bar", ax=ax, fontsize=8)
ax.set_xlabel("Annotation Type")
ax.set_ylabel("# PSM")
plt.savefig("1_annotation_type.pdf", transparent=True)

# plot PTMs per Classification
min_occurence = 10
dup_cols = ['Sequence', 'sample', 'pride id', 'Annotated Mod']
mod_count = df_ann[df_ann["Annotation Type"]
                   == "Fully"].drop_duplicates(dup_cols)
mod_count = mod_count.groupby("Annotated Mod Classification").agg({
    "Annotated Mod": "value_counts"})
mod_count = mod_count.reset_index(level=0)
nrows = mod_count["Annotated Mod Classification"].nunique()
fig, axes = plt.subplots(nrows, 1, figsize=(10, 20), constrained_layout=True)
axes = axes.flatten()
index = 0
for group, ptm_count in mod_count.groupby("Annotated Mod Classification"):
    ax = axes[index]
    try:
        ptm_count[ptm_count["Annotated Mod"] >= min_occurence].plot(
            kind="bar", ax=ax, legend=False)
    except IndexError:
        continue
    ax.set_xlabel("")
    ax.set_title(group)
    ax.set_yscale("log")
    index += 1

while index < nrows:
    axes[index].remove()
    index += 1

plt.savefig("2_ptms_count.pdf", transparent=True)

# plot PTM sites per Classification
dup_cols = ['Sequence', 'sample', 'pride id', 'Annotated Mod']
mod_count = df_ann[df_ann["Annotation Type"]
                   == "Fully"].drop_duplicates(dup_cols)
mod_count = mod_count.groupby("Annotated Mod Classification").agg({
    "Annotated Mod Site": "value_counts"})
mod_count = mod_count.reset_index(level=0)
nrows = mod_count["Annotated Mod Classification"].nunique()
fig, axes = plt.subplots(nrows, 1, figsize=(10, 20), constrained_layout=True)
axes = axes.flatten()
index = 0
for group, site_count in mod_count.groupby("Annotated Mod Classification"):
    ax = axes[index]
    site_count.plot(kind="bar", ax=ax, legend=False)
    ax.set_xlabel("")
    ax.set_title(group)
    ax.set_yscale("log")
    index += 1

while index < nrows:
    axes[index].remove()
    index += 1

plt.savefig("2_ptmsites_count.pdf", transparent=True)


min_occurence = 1000
tokeep = df_ann["Annotated Mod"].value_counts(
)[df_ann["Annotated Mod"].value_counts() >= min_occurence].index.tolist()
mod_count = df_ann[(df_ann["Annotation Type"] == "Fully") &
                   (df_ann["Annotated Mod"].isin(tokeep))]
mod_count = mod_count.groupby("Annotated Mod Classification").agg({
    "Annotated Mod": "value_counts"})
ptm_types = ["Post-translational", "Other",
             "Multiple", "Chemical derivative", "Artefact"]
mod_count = mod_count.reset_index(level=0)

fig, axes = plt.subplots(1, 5, figsize=(3, 2.5), constrained_layout=True, gridspec_kw={
                         'width_ratios': [5, 1, 2, 6, 4]}, sharey=True)
axes = axes.flatten()
for i in range(len(ptm_types)):
    ax = axes[i]
    cond = mod_count["Annotated Mod Classification"] == ptm_types[i]
    mod_count[cond]["Annotated Mod"].plot(kind="bar", ax=ax, legend=False)
    ax.set_ylabel("")
    ax.set_xlabel("")
    ax.set_title(ptm_types[i])
    ax.set_yscale("log")

plt.savefig("2_ptms_global_count.pdf", transparent=True)

mod_count = mod_count.reset_index(level=0)
nrows = mod_count["Annotated Mod Classification"].nunique()
fig, axes = plt.subplots(3, 2, figsize=(5, 7), constrained_layout=True)
axes = axes.flatten()
index = 0
for group, ptm_count in mod_count.groupby("Annotated Mod Classification"):
    ax = axes[index]
    try:
        ptm_count[ptm_count["Annotated Mod"] >= min_occurence].plot(
            kind="bar", ax=ax, legend=False)
    except IndexError:
        continue
    ax.set_xlabel("")
    ax.set_title(group)
    ax.set_yscale("log")
    index += 1

while index < nrows:
    axes[index].remove()
    index += 1

plt.savefig("2_ptms_global_count.pdf", transparent=True)

fig, axes = plt.subplots(2, 1, figsize=(
    9, 9), constrained_layout=True, sharey=True)
ax = axes[0]
per_ptm = df_loc.merge(df_ann)
per_ptm = per_ptm[per_ptm["Annotated Mod"].isin(ptm_count.index)]
per_ptm = per_ptm.groupby("Annotated Mod").apply(
    lambda x: x["Posterior Probability"]).unstack(level=0)
per_ptm = per_ptm[ptm_count.index]
per_ptm.plot(kind="box", ax=ax, showfliers=False)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylim(0.8, 1)
ax.set_xlabel("PTMs Posterior Probability")
ax.set_ylabel("Posterior Probability")

ax = axes[1]
per_saav = df_loc.merge(df_ann)
per_saav = per_saav[per_saav["Annotated Mod"].isin(saav_count.index)]
per_saav = per_saav.groupby("Annotated Mod").apply(
    lambda x: x["Posterior Probability"]).unstack(level=0)
per_saav = per_saav.loc[:, [x in per_saav.columns for x in saav_count.index]]
per_saav.plot(kind="box", ax=ax, showfliers=False)
ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
ax.set_ylim(0.8, 1)
ax.set_xlabel("SAAVs Posterior Probability")
plt.savefig("3_ptm_PER.pdf", transparent=True)

fig, ax = plt.subplots(2, 3, figsize=(9, 6), constrained_layout=True)
ax = ax.flatten()
cond1 = (df_ann["Annotation Type"] == "Fully")
cond2 = (df_ann["Annotated Mod"].isin(ptm_count.index))
position_count = df_ann[cond1 & cond2].drop_duplicates(
    ['Sequence', 'sample', 'pride id', 'Annotated Mod'])

position_count = position_count.groupby(position_count.Sequence.str.len()).apply(
    lambda x: x["Position"].value_counts()).loc[[8, 9, 10, 11, 12], :].unstack(level=1)
position_count = position_count.loc[:, position_count.columns.sort_values()]

for ax_i, i in enumerate([8, 9, 10, 11, 12]):
    plot_pos = position_count.loc[i, 0:i+1]
    plot_pos.index = plot_pos.index.map(lambda x: str(int(x))).str.replace(
        "^0$", "N-term").str.replace(str(i+1), "C-term")
    plot_pos.plot(kind="bar", ax=ax[ax_i])

plt.savefig("4_PTM_position.pdf", transparent=True)


# plot PTMs per sample
with PdfPages('5_studies_PTMs_highcount.pdf') as pdf:
    prev_study = ""
    ncols = 4
    nrows = 3
    max_index = nrows * ncols
    for pride_id, df1 in df_ann.groupby("pride id"):
        fig, axes = gen_figure(nrows, ncols)
        index = 0
        fig.suptitle(pride_id)
        for sample, df2 in df1.groupby("sample"):
            if index >= max_index:
                pdf.savefig(transparent=True)
                plt.close()
                fig, axes = gen_figure(nrows, ncols)
                index = 0
            try:
                df2 = df2.groupby("Dataset Name").agg(
                    {"Annotated Mod": "value_counts"}).unstack("Dataset Name")
            except IndexError:
                print(f"skipping: {sample}")
                continue
            plot_df = df2[df2 > 10].dropna()
            if not plot_df.empty:
                plot_df.plot(kind="bar", ax=axes[index], legend=False)
            else:
                print(f"skipping: {sample}")
                continue
            plt.xticks(rotation=90)
            axes[index].set_title(sample)
            axes[index].set_xlabel("")
            index += 1
        while index < max_index:
            axes[index].remove()
            index += 1
        pdf.savefig(transparent=True)
        plt.close()


with PdfPages('6_studies_PTMs_lowcount.pdf') as pdf:
    prev_study = ""
    ncols = 4
    nrows = 3
    max_index = nrows * ncols
    for pride_id, df1 in df_ann.groupby("pride id"):
        fig, axes = gen_figure(nrows, ncols)
        index = 0
        fig.suptitle(pride_id)
        for sample, df2 in df1.groupby("sample"):
            if index >= max_index:
                pdf.savefig(transparent=True)
                plt.close()
                fig, axes = gen_figure(nrows, ncols)
                index = 0
            try:
                df2 = df2.groupby("Dataset Name").agg(
                    {"Annotated Mod": "value_counts"}).unstack("Dataset Name")
            except IndexError:
                print(f"skipping: {sample}")
                continue
            plot_df = df2[df2 <= 10].dropna()
            if not plot_df.empty:
                plot_df.plot(kind="bar", ax=axes[index], legend=False)
            else:
                print(f"skipping: {sample}")
                continue
            plt.xticks(rotation=90)
            axes[index].set_title(sample)
            axes[index].set_xlabel("")
            index += 1
        while index < max_index:
            axes[index].remove()
            index += 1
        pdf.savefig(transparent=True)
        plt.close()


# unknown mass shifts
df_ann["Mass Shift bin"] = pd.cut(
    df_ann["Mass Shift"], bins=np.arange(-150, 300, 0.02))
partial = df_ann[df_ann["Annotation Type"] == "Partially"].groupby(
    "AA").agg(PSM=("Mass Shift bin", "value_counts"))
partial["Annotation Type"] = "Partially"
unknown = df_ann[df_ann["Annotation Type"] == "None"].groupby(
    "AA").agg(PSM=("Mass Shift bin", "value_counts"))
unknown["Annotation Type"] = "None"

df_unknown = pd.concat([partial[partial.PSM >= 1000],
                        unknown[unknown.PSM >= 1000]])
df_unknown = df_unknown.reset_index().set_index(["Mass Shift bin", "AA"])
df_unknown = df_unknown.sort_index(level=0)
df_unknown.index = df_unknown.index.map(
    lambda x: str(x[0]) + " (" + str(x[1]) + ')')
df_unknown = df_unknown.reset_index()

cond1 = df_unknown[df_unknown["Annotation Type"] == "None"]["PSM"]
cond2 = df_unknown[df_unknown["Annotation Type"] == "Partially"]["PSM"]
fig, ax = plt.subplots(1, 1, figsize=(2, 1.5), constrained_layout=True)
plt1 = ax.barh(cond1.index, cond1, color="blue")
plt2 = ax.barh(cond2.index, cond2, color="orange")
ax.legend([plt1, plt2], ["Partially", "None"])
ax.set_yticks(df_unknown.index)
ax.set_yticklabels(df_unknown["index"])
plt.savefig("unknown_PTMs.pdf", transparent=True)

# circle plot
all_spectra = set(df_fil["Spectrum Name"])
all_mod_spectra = set(df_ann["Spectrum Name"])
fully_spec = set(df_ann[df_ann["Annotation Type"] == "Fully"]["Spectrum Name"])
partially_spec = set(
    df_ann[df_ann["Annotation Type"] == "Partially"]["Spectrum Name"])
unknown_spec = set(df_ann[df_ann["Annotation Type"]
                          == "None"]["Spectrum Name"])
spectra_prop = [len(all_spectra.difference(all_mod_spectra)), len(
    fully_spec), len(partially_spec), len(unknown_spec)]
spectra_labels = ["No PTM", "Fully", "Partially", "None"]


def pct(value, allvals):
    pct = value / np.sum(allvals) * 100
    return " {:.2f}%".format(pct)


for i, label in enumerate(spectra_labels):
    spectra_labels[i] = label + pct(spectra_prop[i], spectra_prop)

# Give color names
fig, ax = plt.subplots(1, 1, figsize=(3, 3), constrained_layout=True)
my_circle = plt.Circle((0, 0), 0.7, color='white')
ax.pie(spectra_prop, labels=spectra_labels, explode=[0, 0.3, 0.3, 0.3])
ax.add_artist(my_circle)
plt.savefig("piechart.pdf", transparent=True)


noptm_aa_count = df_fil[~df_fil["Spectrum Name"].isin(
    df_ann["Spectrum Name"])]["Sequence"].map(list).explode().value_counts()
ptm_aa_count = df_ann["AA"].value_counts()
df_aa_count = pd.concat([noptm_aa_count, ptm_aa_count], axis=1)
df_aa_count = df_aa_count.loc[df_aa_count.sum(
    axis=1).sort_values(ascending=False).index, :]
df_aa_count.columns = ["No PTM", "PTM"]
fig, ax = plt.subplots(1, 1, figsize=(2, 1), constrained_layout=True)
df_aa_count.plot(kind="bar", ax=ax, legend=True,
                 color=["b", "k"], stacked=True)
plt.savefig("os_aa_count.pdf", transparent=True)
