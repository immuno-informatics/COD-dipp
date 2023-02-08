"""
    verify open seach SAAVs if detected by denovo
    fetch SAAVs genomic locations
    add COSMIC and dbSNP annotation to the verified SAAVs
"""
import os
from glob import glob

import numpy as np
import pandas as pd
from Bio import SeqIO


def get_prot_pos(x, db):
    pos = x["Modification peptide position"]
    if pos <= 0:
        pos = 1
    if pos > len(x["Open search peptide"]):
        pos = len(x["Open search peptide"])
    return db[x["Protein Accession"].split(" ")[0]].index(
        x["Open search peptide"]) + pos


def validation_fun(df, db):
    peptide = df["Denovo peptide"]
    peptide_os = df["Open search peptide"]
    splits = df["Annotated Modification"].split(" ")
    position = int(df["Modification peptide position"])
    if position <= 0:
        position = 1
    if position > len(peptide_os):
        position = len(peptide_os)
    try:
        aa1 = aa_conv[splits[0].split("->")[0]]
        aa2 = aa_conv[splits[0].split("->")[1]]
    except KeyError:
        return (False, np.nan)
    if isinstance(peptide, str) and peptide != "":
        if len(peptide) == len(peptide_os):
            validated = peptide[position - 1] == aa2
        else:
            validated = np.nan
    else:
        validated = np.nan
    return (validated, aa1 + str(df["Modification protein position"]) + aa2)


def hamming(s1, s2):
    result = 0
    if len(s1) != len(s2):
        return
    else:
        for x, (i, j) in enumerate(zip(s1, s2)):
            if i != j:
                result += 1
    return result


def read_fasta(fasta_file):
    dict_seq = {}
    for bioseq in SeqIO.parse(fasta_file, format="fasta"):
        dict_seq[bioseq.name] = str(bioseq.seq)
    return dict_seq


def uniprot_var_parse(x, aa_conv):
    try:
        return aa_conv[x[2:5]] + x[5:-3] + aa_conv[x[-3:]]
    except KeyError:
        return np.nan


def isambigious(x):
    if all(x):
        return False
    elif any(x):
        return True
    return np.nan


# dicts
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
    "Leu/Ile": "L",
    "Xle": "L"
}

aa_conv_rev = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "E": "Glu",
    "Q": "Gln",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val"
}

# yapf: disable
# inputs
file_anno = "opensearch_PTMiner_anno_HLA-I.tsv.gz"

# annotation inputs
file_samples = "sample_centered_table.tsv"
# link: 10.6084/m9.figshare.16538097 in pipeline_annotation_files.zip
database_file = "2019-04-30-td-Homo_sapiens_GRCh38_biomart.fasta"
file_denovo = "denovo90ALC_exon_spectra_3m.tsv"  # provided upon request
file_cosmic = "PATH/TO/COSMIC/grch38_v92/CosmicMutantExport.tsv"
file_dbSNP = "PATH/TO/UniProt/homo_sapiens_variation.txt"
# in article supp data

files_dbNSFP = "PATH/TO/dbNSFP/dbNSFP4.1a_variant_proteincentric.chr*"

# outputs
file_output1 = "saavs.tsv"
file_output2 = "summary_global.tsv"
file_output3 = "summary_by_sample.tsv"
file_output4 = "saavs_dbNSFP.tsv"
# yapf: enable

print("reading Open Search and Denovo DataFrames")
dn_cols = [
    'denovo_seq_nomod', 'predicted_score', 'feature_id', 'pride id', 'sample'
]
os_cols = [
    'Spectrum Name', 'Sequence', 'Charge', 'Protein Access', 'Annotated Mass',
    'Annotated Mod', 'Position', 'pride id', 'sample'
]
df_os = pd.read_csv(file_anno, sep="\t", usecols=os_cols)
df_samples = pd.read_csv(file_samples, sep="\t")
df_dn = pd.read_csv(file_denovo, sep="\t", header=0, usecols=dn_cols)

print("reading dbSNP and COSMIC DataFrames")
cosmic_cols = ["GENOMIC_MUTATION_ID", "HGVSP"]
df_cosmic = pd.read_csv(file_cosmic, sep="\t", header=0, usecols=cosmic_cols)
df_cosmic = df_cosmic.drop_duplicates("HGVSP")
df_cosmic.columns = df_cosmic.columns.str.replace("GENOMIC_MUTATION_ID",
                                                  "COSMIC id")
df_cosmic["HGVSP"] = df_cosmic["HGVSP"].str.replace("\.[0-9]+:p\.", ":p.")
df_cosmic = df_cosmic.dropna(
    axis=0, how="any").groupby("HGVSP").agg({
        "COSMIC id": lambda x: ", ".join(x)
    }).reset_index()

dbsnp_cols = [
    "Ensembl translation ID          ", "Variant AA Change   ",
    "Source DB ID    "
]
df_uniprot = pd.read_csv(
    file_dbSNP, sep="\t", skiprows=161, usecols=dbsnp_cols)
df_uniprot = df_uniprot.iloc[1:].drop_duplicates()
col_conv = {
    "Ensembl translation ID          ": "Protein ID",
    "Variant AA Change   ": "Variation",
    "Source DB ID    ": "dbSNP id"
}
df_uniprot.columns = df_uniprot.columns.map(col_conv)
df_uniprot = df_uniprot[df_uniprot["Protein ID"].str.contains(
    "^ENSP", na=False)]
df_uniprot["Variation"] = df_uniprot["Variation"].map(
    lambda x: uniprot_var_parse(x, aa_conv))
df_uniprot = df_uniprot[~df_uniprot["Variation"].isnull()].drop_duplicates()
df_uniprot = df_uniprot.groupby(["Variation", "Protein ID"]).agg({
    "dbSNP id":
    lambda x: ", ".join(x)
}).reset_index()

print("filtering dataframes")
# excluding HLA II samples
HLAI_samples = df_samples[df_samples["HLA Class"] == 1]["Sample Name"].tolist()
df_os = df_os[df_os["sample"].isin(HLAI_samples)].reset_index(drop=True).copy(
    deep=True)

# filter dataframes

col_conv_os = {
    'Spectrum Name': 'Spectrum Name',
    'Protein Access': 'Protein Accession',
    'Sequence': 'Open search peptide',
    'Annotated Mod': 'Annotated Modification',
    'Annotated Mass': 'Delta mass',
    'pride id': 'Pride ID',
    'sample': 'Sample ID',
    'Position': 'Modification peptide position'
}

col_conv_dn = {
    'denovo_seq_nomod': 'Denovo peptide',
    'predicted_score': 'Denovo score',
    'spectrum': 'Spectrum Name',
    'Annotated Mod': 'Annotated Modification',
    'Annotated Mass': 'Delta mass',
    'pride id': 'Pride ID',
    'sample': 'Sample ID',
    'Charge': 'Charge'
}

df_dn["spectrum"] = df_dn["feature_id"].map(lambda x: x.rsplit(":", 1)[0])
df_dn.drop("feature_id", axis=1, inplace=True)
df_os.columns = df_os.columns.map(col_conv_os)
df_dn.columns = df_dn.columns.map(col_conv_dn)

df_os["SAAV"] = df_os["Annotated Modification"].str.contains(
    "->", regex=False, na=False)

saav_specs = df_os[~df_os["Annotated Modification"].isnull()].groupby(
    "Spectrum Name").agg({
        "SAAV": isambigious
    })
df_os = df_os.drop("SAAV", axis=1)
saav_specs.reset_index(inplace=True)
saav_specs.columns = ["Spectrum Name", "Ambiguous mass shift"]
saav_specs = saav_specs[~saav_specs["Ambiguous mass shift"].isnull()]

df_os = df_os[df_os["Annotated Modification"].str.contains("->", na=False)]
df_os = df_os.merge(saav_specs, how="left")

df_os = df_os[~df_os["Protein Accession"].str.contains("^sp")].reset_index(
    drop=True).copy(deep=True)

# read the protein database to retrieve the peptides
# location in the proteins
print("retreiving SAAVs location within proteins")

db = read_fasta(database)
df_os["Modification protein position"] = df_os.apply(
    lambda x: get_prot_pos(x, db), axis=1)

print("verifying if SAAVs detected by denovo")
df = df_dn.merge(df_os, how="right")

print("verifying SAAVs with denovo")
# replace I by L since they have the same isotopic mass
# denovo cannot discriminate
os_peptides = df["Open search peptide"].str.replace("I", "L")
dn_peptides = df["Denovo peptide"].fillna("")
dn_peptides = dn_peptides.str.replace("I", "L")

temp_df = pd.concat([os_peptides, dn_peptides], axis=1)

df["Mismatch with denovo"] = temp_df.apply(
    lambda x: hamming(x["Open search peptide"], x["Denovo peptide"]
                      ) if x["Denovo peptide"] != "" else np.nan,
    axis=1)
del temp_df

df_saav = df.reset_index(drop=True).copy(deep=True)
df_saav["Modification peptide position"] = df_saav[
    "Modification peptide position"].map(int)
df_saav["Modification protein position"] = df_saav[
    "Modification protein position"].map(int)
values = df_saav.apply(
    lambda x: validation_fun(x, db), axis=1).apply(pd.Series)
df_saav["Variation"] = values[1]
df_saav["Validated with denovo"] = values[0]

df_saav = df_saav[~df_saav.Variation.isnull()].reset_index(drop=True).copy(
    deep=True)
df_saav["Gene Name"] = df_saav["Protein Accession"].map(
    lambda x: x.split(" ", 1)[0].split("|")[3])
df_saav["Gene ID"] = df_saav["Protein Accession"].str.split(
    "|", expand=True)[2]
df_saav["Transcript ID"] = df_saav["Protein Accession"].str.split(
    "|", expand=True)[1]
df_saav["Protein ID"] = df_saav["Protein Accession"].str.split(
    "|", expand=True)[0]
df_saav["HGVSP"] = df_saav["Protein ID"] + ":p." + df_saav.Variation.map(
    lambda x: aa_conv_rev[x[0]] + x[1:-1] + aa_conv_rev[x[-1]])

df_saav = df_saav.merge(
    df_uniprot,
    how="left",
    left_on=["Protein ID", "Variation"],
    right_on=["Protein ID", "Variation"])
df_saav = df_saav.merge(
    df_cosmic, how="left", left_on="HGVSP", right_on="HGVSP")

# generate a summary table
print("generating summary dataframe")
df_sum_global = df_saav[[
    "Gene Name", "Gene ID", "Transcript ID", "Protein ID", "HGVSP",
    "Variation", "Ambiguous mass shift"
]].value_counts().reset_index()
cols = df_sum_global.columns.tolist()
cols[-1] = "#PSM"
df_sum_global.columns = cols

sum_cols1 = [
    "Gene Name", "Gene ID", "Transcript ID", "Protein ID", "HGVSP",
    "Variation", "Ambiguous mass shift", "Sample ID", "Pride ID"
]
sum_cols2 = [
    "Gene Name", "Gene ID", "Transcript ID", "Protein ID", "HGVSP",
    "Variation", "Ambiguous mass shift", "#PSM", "Sample ID", "Pride ID"
]
df_sum_sample = df_saav.groupby(sum_cols1).agg({
    "Spectrum Name": "count"
}).reset_index()
df_sum_sample.columns = df_sum_sample.columns.str.replace(
    "Spectrum Name", "#PSM")
df_sum_sample = df_sum_sample[sum_cols2]
df_sum_sample = df_sum_sample.sort_values(["Sample ID", "Protein ID"])

df_sum_global = df_sum_global.merge(
    df_uniprot,
    how="left",
    left_on=["Protein ID", "Variation"],
    right_on=["Protein ID", "Variation"])
df_sum_global = df_sum_global.merge(
    df_cosmic, how="left", left_on="HGVSP", right_on="HGVSP")

df_sum_sample = df_sum_sample.merge(
    df_uniprot,
    how="left",
    left_on=["Protein ID", "Variation"],
    right_on=["Protein ID", "Variation"])
df_sum_sample = df_sum_sample.merge(
    df_cosmic, how="left", left_on="HGVSP", right_on="HGVSP")

saav_cols = [
    'Spectrum Name', 'Denovo peptide', 'Denovo score', 'Open search peptide',
    'Mismatch with denovo', 'Validated with denovo', 'Protein Accession',
    'Delta mass', 'Annotated Modification', 'Modification peptide position',
    'Modification protein position', 'Variation', 'Gene Name', 'Gene ID',
    'Transcript ID', 'Protein ID', 'HGVSP', 'Ambiguous mass shift', 'dbSNP id',
    'COSMIC id', 'Pride ID', 'Sample ID'
]

df_saav = df_saav[saav_cols]

# dbNSFP annotation
print("Adding dbNSFP annotation")
files_dbNSFP = glob(files_dbNSFP)

df_list = []
for file_ in files_dbNSFP:
    df_dbNSFP = pd.read_csv(file_, sep="\t")
    df_dbNSFP["HGVSP"] = df_dbNSFP["Ensembl_proteinid"] + ":p." + df_dbNSFP[
        "aaref"].map(
            lambda x: aa_conv_rev.get(x, "")) + df_dbNSFP["aapos"].map(
                str) + df_dbNSFP["aaalt"].map(lambda x: aa_conv_rev.get(x, ""))
    df_dbNSFP = df_dbNSFP.drop(
        ["Ensembl_proteinid", "aaref", "aapos", "aaalt"], axis=1)
    df_list.append(df_saav.merge(df_dbNSFP, how="inner"))

df_saav_dbNSFP = pd.concat(df_list)

if not df_saav_dbNSFP.empty:
    hvsp_dbNSFP = df_saav_dbNSFP.set_index("HGVSP").apply(
        lambda x: sum(x[x.index.str.contains("rankscore")] >= 0.75),
        axis=1).reset_index().drop_duplicates()
    hvsp_dbNSFP.columns = ["HGVSP", "dbNSFP rankscore >= 0.75"]
    hvsp_dbNSFP = hvsp_dbNSFP.sort_values(
        "dbNSFP rankscore >= 0.75", ascending=False).drop_duplicates("HGVSP")
    df_saav = df_saav.merge(hvsp_dbNSFP, how="left")
    df_sum_global = df_sum_global.merge(hvsp_dbNSFP, how="left")
    df_sum_sample = df_sum_sample.merge(hvsp_dbNSFP, how="left")

print("writing output")
df_saav.to_csv(file_output1, sep="\t", header=True, index=False)
df_sum_global.to_csv(file_output2, sep="\t", header=True, index=False)
df_sum_sample.to_csv(file_output3, sep="\t", header=True, index=False)
df_saav_dbNSFP.to_csv(file_output4, sep="\t", header=True, index=False)

# t = df_saav.groupby(df_saav.Variation.map(lambda x: x[0] + " -> " + x[-1])).apply(lambda x: any(x["Ambiguous mass shift"]))
# t[t]
# (A -> E) can be confused with other PTMs
# (A -> N) can be confused with other PTMs
# (A -> Q) can be confused with other PTMs
# (A -> V) can be confused with other PTMs
# (D -> E) can be confused with other PTMs
# (D -> W) can be confused with other PTMs
# (E -> W) can be confused with other PTMs
# (F -> Y) can be confused with other PTMs
# (G -> N) can be confused with other PTMs
# (G -> Q) can be confused with other PTMs
# (K -> R) can be confused with other PTMs
# (L -> R) can be confused with other PTMs
# (M -> F) can be confused with other PTMs
# (N -> D) can be confused with other PTMs
# (N -> Q) can be confused with other PTMs
# (P -> E) can be confused with other PTMs
# (Q -> E) can be confused with other PTMs
# (S -> A) can be confused with other PTMs
# (S -> D) can be confused with other PTMs
# (S -> E) can be confused with other PTMs
# (S -> T) can be confused with other PTMs
# (T -> E) can be confused with other PTMs
# (V -> L) can be confused with other PTMs
# (V -> R) can be confused with other PTMs
# (V -> W) can be confused with other PTMs
