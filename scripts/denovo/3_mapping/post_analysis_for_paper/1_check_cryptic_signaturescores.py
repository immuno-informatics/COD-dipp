"""
    Checks if Denovo cryptic peptides singature scores using the closed search PSWMs
"""

import os
import pickle as pkl
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

import numpy as np


def runcommand(cmd):
    proc = subprocess.Popen(cmd,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True,
                            universal_newlines=True)
    std_out, std_err = proc.communicate()
    return proc.returncode, std_out, std_err


# yapf: disable
root = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all"
file_psmws = os.path.join(root, "signatures_analysis/reports/pswm_umap.pkl")
file_peptides = os.path.join(root, "peptides_probamconvert/gene_inference/bam_files_2/coverage/merged.sorted_features.tsv")

scoring_script = "$HOME/git_repositories/ipip/scripts/denovo/3_mapping/post_analysis_for_paper/pswm2score.py"

root_output = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation"
output_exo_scores = os.path.join(root_output, "overlap/exonic_PSWMscores.tsv")
output_exo_bestscores = os.path.join(root_output, "overlap/exonic_bestPSMWscores.tsv")
output_cryp_scores = os.path.join(root_output, "overlap/cryptic_PSWMscores.tsv")
output_cryp_bestscores = os.path.join(root_output, "overlap/cryptic_bestPSWMscores.tsv")
pdf_output = os.path.join(root_output, "overlap/cryptic_pswmScores.pdf")
# yapf: enable

# loading data
df_all = pd.read_csv(file_peptides, sep="\t")
df = df_all[df_all["tag"].str.contains("denovo")].copy(deep=True)
df.columns = df.columns.str.replace("sample_name", "Sample Name")
df.drop("study_id", axis=1, inplace=True)

fasta_input = "/tmp/temp_cryptic.fasta"
output = "/tmp/temp_cryptic_scores.tsv"
dfscores_list = []
cond = df["tag"] == "denovo cryptic"
for sample_name, sdf in df[cond].groupby("Sample Name"):
    print(f"processing intronic {sample_name}")
    command = f"python {scoring_script} "
    command += f"-input_fasta {fasta_input} -input_pswm {file_psmws} "
    command += f"-output {output} -sample_name '{sample_name}'"
    peps = sdf["peptide"].drop_duplicates().tolist()
    with open(fasta_input, "w") as out:
        for i, pep in enumerate(peps):
            out.write(f">seq_{i}\n{pep}\n")
    code, out, err = runcommand(command)
    if code != 0:
        print(f"Error: {sample_name} scoring failed")
        print(err)
        continue
    df_temp = pd.read_csv(output, sep="\t", header=0)
    df_temp["Sample Name"] = sample_name
    dfscores_list.append(df_temp)

os.system(f"rm {fasta_input} {output}")
df_cryp_scores = pd.concat(dfscores_list, ignore_index=True)

df_cryp_bestscores = df_cryp_scores.groupby("Sample Name").apply(
    lambda x: x.groupby("peptide").apply(lambda z: z.sort_values(
        "normalized score").iloc[-1])).reset_index(drop=True)

dfscores_list = []
cond = df["tag"] != "denovo cryptic"
for sample_name, sdf in df[cond].groupby("Sample Name"):
    print(f"processing exonic sequences for {sample_name}")
    command = f"python {scoring_script} "
    command += f"-input_fasta {fasta_input} -input_pswm {file_psmws} "
    command += f"-output {output} -sample_name '{sample_name}'"
    peps = sdf["peptide"].drop_duplicates().tolist()
    with open(fasta_input, "w") as out:
        for i, pep in enumerate(peps):
            out.write(f">seq_{i}\n{pep}\n")
    code, out, err = runcommand(command)
    if code != 0:
        print(f"Error: {sample_name} scoring failed")
        print(err)
        continue
    df_temp = pd.read_csv(output, sep="\t", header=0)
    df_temp["Sample Name"] = sample_name
    dfscores_list.append(df_temp)

df_exo_scores = pd.concat(dfscores_list, ignore_index=True)

df_exo_bestscores = df_exo_scores.groupby("Sample Name").apply(
    lambda x: x.groupby("peptide").apply(lambda z: z.sort_values(
        "normalized score").iloc[-1])).reset_index(drop=True)

df_exo_scores.to_csv(output_exo_scores, sep="\t", header=True, index=False)
df_exo_bestscores.to_csv(output_exo_bestscores,
                         sep="\t",
                         header=True,
                         index=False)

df_cryp_scores.to_csv(output_cryp_scores, sep="\t", header=True, index=False)
df_cryp_bestscores.to_csv(output_cryp_bestscores,
                          sep="\t",
                          header=True,
                          index=False)

fig, ax = plt.subplots(1, 1, figsize=(4, 3))
data1 = df_exo_bestscores["normalized score"].tolist()
data2 = df_cryp_bestscores["normalized score"].tolist()
values2, bins2, plt2 = ax.hist(
    data2,
    bins=np.arange(-1, 1.1, 0.1),
    alpha=0.75,
    color="red",
    density=True,
    edgecolor='white',
)
values1, bins1, plt1 = ax.hist(data1,
                               bins=np.arange(-1, 1.1, 0.1),
                               edgecolor='white',
                               alpha=0.75,
                               density=True)

handle1 = mpl.lines.Line2D([], [])
handle2 = mpl.lines.Line2D([], [], c='red')
ax.legend([handle1, handle2],
          ["Exonic denovo peptides", "Cryptic denovo peptides"])
shared_area = sum(
    np.diff(bins1) * np.min(np.vstack([values2, values1]), axis=0))
ax.text(-0.75, 1.5, f"shared area: {round(shared_area, 2)}", fontsize=8)
fig.suptitle("Denovo peptides best PSWMs scores", fontsize=8)
plt.savefig(pdf_output, transparent=True)
