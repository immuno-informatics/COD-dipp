import pandas as pd
import os
from sklearn.metrics import r2_score

root= "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/peptides_probamconvert/denovo/output"
file_bindings = os.path.join(root, "overlap/bestcryptic_peptides_bindings.tsv")
file_cryp_rt = os.path.join(root, "overlap/RT_cryptic_testing_data.tsv")
file_exo_rt = os.path.join(root, "overlap/RT_exonic_testing_data.tsv")

df_cryp_rt = pd.read_csv(file_cryp_rt, sep="\t", header=0)
df_bindings = pd.read_csv(file_bindings, sep="\t", header=0)

peptides = df_bindings.Peptide.tolist()
pep_scores = df_cryp_rt.groupby("peptide").agg({"predicted score": "mean"})
pep_thresholds  = df_cryp_rt.groupby("peptide").agg({"RT threshold": "mean"})

pep_stats = pd.concat([pep_thresholds.loc[peptides], pep_scores.loc[peptides]], axis=1)
pep_stats["pass threshold"] = pep_stats["predicted score"] >= pep_stats["RT threshold"]

pep_stats.reset_index(inplace=True)
pep_stats.columns = ['Peptide', 'mean RT threshold (>=0.7 r-squared)', 'mean denovo score', 'pass mean RT threshold']
df_bindings = df_bindings.merge(pep_stats)

df_bindings = df_bindings.merge(pep_stats)
df_bindings.to_csv(file_bindings, sep="\t", header=True, index=False)
