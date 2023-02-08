import numpy as np
import pandas as pd

# change to output of "./1_generate_annotation.py"
tsv_file2 = "CHANGE2FilePath"

df_features = pd.read_csv(tsv_file2, sep="\t")
print(df_features.shape)
print(df_features.head)
global remainder


def extract_frame2(x, ascending):
    global remainder
    if ascending:
        if x["Intron end"] == -1:
            return np.nan
        positions = list(range(x.start - remainder, x.end + 1, 3))
        ex_in_pos = list(range(x.start - remainder, x["Intron end"] + 1, 3))
        if remainder != 0:
            positions = positions[1:]
            ex_in_pos = ex_in_pos[1:]
        try:
            remainder = x.end - positions[-1]
        except IndexError:
            remainder = x.end
    else:
        if x["Intron start"] == -1:
            return np.nan
        positions = list(range(x.end + remainder, x.start - 1, -3))
        ex_in_pos = list(range(x.end + remainder, x["Intron start"] - 1, -3))
        if remainder != 0:
            positions = positions[1:]
            ex_in_pos = ex_in_pos[1:]
        try:
            remainder = positions[-1] - x.start
        except IndexError:
            remainder = x.start
    try:
        return [ex_in_pos[0], ex_in_pos[-1]]
    except IndexError:
        return np.nan


def extract_frame(df):
    if isinstance(df, pd.Series):
        # 1 entry per transcript: case of transcript consisting of 1 exon only
        # i.e: https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=ENSG00000284546;r=11:4329865-4330449;t=ENST00000640302
        return pd.DataFrame()
    global remainder
    remainder = 0
    ascending = df.strand.iloc[0] == "+"
    # print(ascending)
    coord_col = "start" if ascending else "end"
    # keep everything after translation start site
    df = df.sort_values("start", ascending=ascending, ignore_index=True)
    index = df[df.feature == "Exon"].first_valid_index()
    df = df.iloc[index:].reset_index(drop=True)
    df_exons = df[df.feature == "Exon"].reset_index(drop=True).copy(deep=True)
    if coord_col == "start":
        df_exons["Intron end"] = df[df.feature ==
                                    "Intron"].end.reset_index(drop=True)
        df_exons["Intron end"] = df_exons["Intron end"].fillna(-1).astype(int)
    else:
        df_exons["Intron start"] = df[df.feature ==
                                      "Intron"].start.reset_index(drop=True)
        df_exons["Intron start"] = df_exons["Intron start"].fillna(
            -1).astype(int)
    df_exons["inframe positions"] = df_exons.apply(
        lambda x: extract_frame2(x, ascending), axis=1)
    return df_exons


df_features_inframes = df_features.groupby("name").apply(extract_frame)

df_features_inframes["inframe positions start"] = df_features_inframes["inframe positions"].map(
    lambda x: x[0] if isinstance(x, list) else np.nan)
df_features_inframes["inframe positions end"] = df_features_inframes["inframe positions"].map(
    lambda x: x[1] if isinstance(x, list) else np.nan)
df_features_inframes[["name", "chrom", "start", "end", "feature", "Intron start", "Intron end", "inframe positions start",
                      "inframe positions end"]].to_csv("df_features_inframes.tsv", sep="\t", header=True, index=False)
