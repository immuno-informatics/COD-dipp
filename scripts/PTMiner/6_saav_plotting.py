import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file_input1 = "saavs.tsv"
file_input2 = "summary_global.tsv"
file_input3 = "summary_by_sample.tsv"
file_input4 = "sample_centered_table.tsv"  # in article supp data

df_saav = pd.read_csv(file_input1, sep="\t")
df_sum_global = pd.read_csv(file_input2, sep="\t")
df_sum_sample = pd.read_csv(file_input3, sep="\t")
df_ann = pd.read_csv(file_input4, sep="\t")

# dataframes pre-processing
df_saav["Validated with denovo"].fillna(False, inplace=True)
df_ann.columns = df_ann.columns.str.replace("Sample Name", "Sample ID")
df_saav = df_saav.merge(df_ann, how="left")
df_sum_sample = df_sum_sample.merge(df_ann, how="left")

# plotting
# barplot of the validated SAAVs count
i = 1
fig, ax = plt.subplots(1, 1, figsize=(1.5, 2), constrained_layout=True)
matched_denovo = df_saav.groupby(["Protein ID", "Variation"]).apply(
    lambda x: any(~x["Denovo score"].isnull())).value_counts()
non_matched_denovo = matched_denovo[False]
matched_denovo = matched_denovo[True]
validated_denovo = df_saav.groupby(["Protein ID", "Variation"]).apply(
    lambda x: any(x["Validated with denovo"])).value_counts()[True]
non_validated_denovo = matched_denovo - validated_denovo
df_count = pd.DataFrame({"Validated": [validated_denovo, np.nan], "Non validated": [
                        non_validated_denovo, non_matched_denovo]}, index=["Matched", "Non matched"])

df_count.plot(kind="bar", legend=True, ax=ax, stacked=True)
for x, y in enumerate(df_count.sum(axis=1)):
    ax.text(x, y + 10, str(int(y)), ha="center", fontsize=8)

ax.set_ylabel("# of SAAVs")
fig.suptitle("denovo agreement", fontsize=8)
plt.savefig("plot{}.pdf".format(i), transparent=True)

# number of SAAVs per sample
i = 2
fig, ax = plt.subplots(1, 1, figsize=(2, 2), constrained_layout=True)
df_count_per_sample = df_saav.groupby(["Sample ID", "Pride ID"]).apply(
    lambda x: x[["Protein ID", "Variation"]].drop_duplicates().shape[0])
df_count_per_sample.plot(kind="box", legend=False, ax=ax, showfliers=False)
ax.set_xlabel("# SAAVs per sample")
ax.set_xticklabels([])
plt.savefig("plot{}.pdf".format(i), transparent=True)

# number of SAAVs per condition
i = 3
fig, ax = plt.subplots(1, 1, figsize=(4, 1), constrained_layout=True)
agg_funs = {"COSMIC id": lambda z: any(
    ~z.isnull()), "dbSNP id": lambda z: any(~z.isnull())}
df_count_per_condition = df_sum_sample.groupby("Disease").apply(
    lambda x: x.groupby(["Protein ID", "Variation"]).agg(agg_funs))
df_count_per_condition["dbSNP id"] = df_count_per_condition.apply(
    lambda x: (x["dbSNP id"] == True) & (x["COSMIC id"] == False), axis=1)
df_count_per_condition = df_count_per_condition.reset_index()

df_count_per_condition["Unreported SAAV"] = df_count_per_condition.apply(
    lambda x: x["COSMIC id"] == False and x["dbSNP id"] == False, axis=1)

df_count_per_condition = df_count_per_condition.groupby("Disease").agg(sum)
df_count_per_condition.columns = ['COSMIC', 'dbSNP', 'Unreported SAAV']
index = df_count_per_condition.sum(axis=1).sort_values().index.tolist()
index.remove("Na")
# nsamples = df_ann.groupby("Disease")["Sample Name"].agg("count").sort_values()
# df_count_per_condition = df_count_per_condition.apply(lambda x: x/nsamples[x.name], axis=1)
index_to_plot = [index[-1], index[-3], index[-5],
                 index[-6], index[-7], index[-8], index[-10]]
df_count_per_condition.loc[index_to_plot, :].plot(
    kind="barh", ax=ax, stacked=True, legend=True)
ax.set_xscale("log")
plt.savefig("plot{}.pdf".format(i), transparent=True)

# plot SNV vs NON SNV (DbSNP, COSMIC, Gene Census)
i = 4
fig, ax = plt.subplots(1, 1, figsize=(2, 3), constrained_layout=True)
df_sum_global["dbSNP id"] = df_sum_global.apply(
    lambda x: np.nan if isinstance(x["COSMIC id"], str) else x["dbSNP id"], axis=1)
mutation_type = df_sum_global.groupby("Ambiguous mass shift").agg(
    {"dbSNP id": "count", "COSMIC id": "count", "HGVSP": "count"}).T
mutation_type.index = mutation_type.index.str.replace(
    "HGVSP", "Unreported\nSAAV")

mutation_type.plot(kind="bar", legend=True, ax=ax,
                   stacked=True, color=["#9a5ea1", "#98823c"])
ax.set_xlabel("Mutation type")
ax.set_ylabel("Mutation Log count")
plt.savefig("plot{}.pdf".format(i), transparent=True)
