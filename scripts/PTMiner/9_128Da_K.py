import pandas as pd
import seaborn as sb
import matplotlib.pyplot as plt

file_in = "opensearch_PTMiner_loc_HLA-I.tsv"

df = pd.read_csv(file_in, sep="\t", header=0)

# isolate -128Da on K
df_f = df[(df["Mass shift"] >= -128.1) & (df["Mass shift"] <= -128.08) &
          (df["AA"] == "K")]

# peptide length distribution
df_plot1 = df_f.Sequence.str.len().value_counts().sort_index()

# amino acid position of the mass shift distribution
df_plot2 = df_f.groupby(df_f.Sequence.str.len())["Position"].value_counts(
    normalize=True).unstack(level=1).fillna(0)

# keep same order fr df_plot1 and df_plot2
df_plot2 = df_plot2.loc[df_plot1.index[::-1]]

fig, axes = plt.subplots(1,
                         2,
                         figsize=(6, 3),
                         constrained_layout=True,
                         gridspec_kw={'width_ratios': [1, 10]})

# barplot of peptide length
i = 0
ax = axes[i]
df_plot1.plot(kind="barh", ax=ax)
ax.set_ylabel("Peptide length")

# heatmap of mass shift localization
i += 1
ax = axes[i]
sb.heatmap(df_plot2, ax=ax, cmap="YlGnBu")
ax.set_xlabel("Mass shift amino acid position")
plt.minorticks_off()
ax.set_yticklabels([])
plt.savefig("K_128Da.pdf", transparent=True)
