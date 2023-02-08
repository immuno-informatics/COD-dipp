import pandas as pd
from pyteomics import achrom
import pyteomics
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import r2_score
import pickle as pkl
import os
from matplotlib.backends.backend_pdf import PdfPages


def achrom_exception(x, RCs):
    try:
        return achrom.calculate_RT(x, RCs)
    except pyteomics.auxiliary.structures.PyteomicsError:
        return 0


def predict(df, RCs):
    df_list = []
    for fraction, sdf in df.groupby("fraction"):
        try:
            sdf["RT pred"] = sdf.peptide.map(
                lambda x: achrom_exception(x, RCs[fraction]))
            df_list.append(sdf)
        except KeyError:
            continue
    return pd.concat(df_list, ignore_index=True)


def predict_per_sample(df_cryp, df_exo, df_cs, sample_name):
    df_cs = df_cs[df_cs["Sample Name"] == sample_name].copy(deep=True)
    df_cryp = df_cryp[df_cryp["Sample Name"] == sample_name].copy(deep=True)
    df_exo = df_exo[df_exo["Sample Name"] == sample_name].copy(deep=True)
    if df_cs.empty or df_cryp.empty or df_exo.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(
        ), {}
    df_cryp["fraction"] = df_cryp["spectrum"].map(lambda x: x.split(".")[0])
    df_exo["fraction"] = df_exo["spectrum"].map(lambda x: x.split(".")[0])
    df_training = df_cs.sample(int(df_cs.shape[0] * 0.8))
    df_testing = df_cs[~df_cs.index.isin(df_training.index)]
    df_testing = df_testing[~df_testing.peptide.isin(df_training.peptide
                                                     )].reset_index().copy(
                                                         deep=True)
    RCs = {}
    for fraction, sdf in df_training.groupby("fraction"):
        RCs[fraction] = achrom.get_RCs(sdf.peptide.tolist(),
                                       sdf["RT exp"].tolist(),
                                       lcp=0.3)
    df_testing = predict(df_testing, RCs)
    df_cryp = predict(df_cryp, RCs)
    df_exo = predict(df_exo, RCs)
    return df_cryp, df_exo, df_testing, df_training, RCs


def estimate_threshold(df):
    score_threshold = 0
    r2_squared = r2_score(df["RT exp"], df["RT pred"])
    if r2_squared >= 0.7:
        return score_threshold
    else:
        for i in range(0, 101):
            sdf = df[df["predicted score"] >= i]
            try:
                if sdf.shape[0] < 2:
                    continue
                r2_squared = r2_score(sdf["RT exp"], sdf["RT pred"])
                if r2_squared >= 0.7:
                    score_threshold = i
                    break
            except ValueError:
                continue
        if score_threshold == 0:
            score_threshold = -1
    return score_threshold


def plot(df_testing, df_exo, df_cryp, score_threshold=True):
    if score_threshold:
        exo_score_threshold = estimate_threshold(df_exo)
        cryp_score_threshold = estimate_threshold(df_cryp)
    else:
        exo_score_threshold = 0
        cryp_score_threshold = 0
    fig, axes = plt.subplots(1,
                             3,
                             figsize=(8, 3),
                             constrained_layout=True,
                             sharey=True)
    ymax = max(df_testing["RT exp"].max(), df_exo["RT exp"].max(),
               df_cryp["RT exp"].max())
    axes[0].set_ylim(0, ymax)
    ax = axes[0]
    df_testing.plot(kind="scatter",
                    x="RT exp",
                    y="RT pred",
                    s=0.1,
                    ax=ax,
                    rasterized=True)
    ax.plot(np.linspace(0, 100, 101), np.linspace(0, 100, 101), "k-")
    ax.set_title("closed search RT exp VS pred", fontsize=8)
    ax.annotate("r-squared = {:.3f}".format(
        r2_score(df_testing["RT exp"], df_testing["RT pred"])),
                (0, ax.get_ylim()[1] * 0.95),
                fontsize=8)
    ax = axes[1]
    df_exo2 = df_exo[df_exo["predicted score"] >= exo_score_threshold]
    df_exo2.plot(kind="scatter",
                 x="RT exp",
                 y="RT pred",
                 s=0.1,
                 ax=ax,
                 rasterized=True)
    ax.plot(np.linspace(0, 100, 101), np.linspace(0, 100, 101), "k-")
    ax.set_title("Denovo exonic peptides\nRT exp VS pred", fontsize=8)
    r2_value = r2_score(df_exo2["RT exp"], df_exo2["RT pred"])
    msg = "r-squared = {:.3f}".format(r2_value)
    if exo_score_threshold > 0:
        msg = msg + f"\nscore threshold = {exo_score_threshold}"
        msg = msg + "\nretained fraction {:.1f}%".format(
            (df_exo2.shape[0] / df_exo.shape[0]) * 100)
    ax.annotate(msg, (0, ax.get_ylim()[1] * 0.95),
                fontsize=8,
                verticalalignment="top")
    ax = axes[2]
    df_cryp2 = df_cryp[df_cryp["predicted score"] >= cryp_score_threshold]
    df_cryp2.plot(kind="scatter",
                  x="RT exp",
                  y="RT pred",
                  s=0.1,
                  ax=ax,
                  rasterized=True)
    ax.plot(np.linspace(0, 100, 101), np.linspace(0, 100, 101), "k-")
    ax.set_title("Denovo cryptic peptides\nRT exp VS pred", fontsize=8)
    r2_value = r2_score(df_cryp2["RT exp"], df_cryp2["RT pred"])
    msg = "r-squared = {:.3f}".format(r2_value)
    if cryp_score_threshold > 0:
        msg = msg + f"\nscore threshold = {cryp_score_threshold}"
        msg = msg + "\nretained fraction {:.1f}%".format(
            (df_cryp2.shape[0] / df_cryp.shape[0]) * 100)
    ax.annotate(msg, (0, ax.get_ylim()[1] * 0.95),
                fontsize=8,
                verticalalignment="top")
    return fig, axes, exo_score_threshold, cryp_score_threshold


def main():
    # yapf: disable
    root = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all"
    file_input_cs = os.path.join(root, "tables/closed_search_PSMs.tsv")
    file_input_dn = os.path.join(root, "tables/denovo_90ALC.tsv")
    file_rt = os.path.join(root, "tables/all_specra_rt_exp.tsv")
    file_input_exo = os.path.join(root, "denovo_annotation/overlap/exonic_bestPSMWscores.tsv")
    file_input_cryp = os.path.join(root, "denovo_annotation/overlap/cryptic_bestPSWMscores.tsv")
    file_ann = os.path.join(root, "tables/sample_centered_table.tsv")

    root_output = os.path.join(root, "denovo_annotation/overlap")
    output_pdf = os.path.join(root_output, "RT_quality_checks.pdf")
    df_training_output = os.path.join(root_output, "RT_training_data.tsv")
    df_testing_output = os.path.join(root_output, "RT_testing_data.tsv")
    df_cryp_output = os.path.join(root_output, "RT_cryptic_testing_data.tsv")
    df_exo_output = os.path.join(root_output, "RT_exonic_testing_data.tsv")
    RCs_output = os.path.join(root_output, "RT_RCs_model.pkl")
    # yapf: enable
    print("Reading closed search data")
    df_cs = pd.read_csv(
        file_input_cs,
        sep="\t",
        header=0,
        usecols=["peptide", "spectrum", "RT exp", "pride id", "sample"])
    print("Reading denovo data")
    df_dn = pd.read_csv(file_input_dn,
                        sep="\t",
                        header=0,
                        usecols=[
                            "feature_id", "denovo_seq_nomod",
                            "predicted_score", "pride id", "sample"
                        ])
    print("Reading cryptic denovo peptides data")
    df_cryp = pd.read_csv(file_input_cryp, sep="\t", header=0)
    print("Reading exonic denovo peptides data")
    df_exo = pd.read_csv(file_input_exo, sep="\t", header=0)
    print("Reading experimental MS retention time values")
    df_rts = pd.read_csv(file_rt, sep="\t", header=0)
    df_ann = pd.read_csv(file_ann, sep="\t")
    df_ann = df_ann[df_ann["HLA Class"] == 1]

    df_dn = df_dn[df_dn["sample"].isin(df_ann["Sample Name"])]
    df_cs = df_cs[df_cs["sample"].isin(df_ann["Sample Name"])]

    df_cs["fraction"] = df_cs["spectrum"].map(lambda x: x.split(".")[0])
    df_dn["spectrum"] = df_dn["feature_id"].map(lambda x: x.rsplit(":", 1)[0])
    df_dn["fraction"] = df_dn["spectrum"].map(lambda x: x.split(".")[0])

    df_dn = df_dn.merge(df_rts, how="inner")
    df_dn.columns = df_dn.columns.str.replace("denovo_seq_nomod", "peptide")
    df_dn.columns = df_dn.columns.str.replace("predicted_score",
                                              "predicted score")
    df_dn.columns = df_dn.columns.str.replace("sample", "Sample Name")
    df_dn["predicted score"] = np.exp(df_dn["predicted score"]) * 100

    df_cs.columns = df_cs.columns.str.replace("sample", "Sample Name")

    df_cryp = df_cryp.merge(df_dn, how="inner")
    df_exo = df_exo.merge(df_dn, how="inner")

    df_cryp["RT exp"] = df_cryp["RT exp"] / 60
    df_exo["RT exp"] = df_exo["RT exp"] / 60

    df_cryp_list = []
    df_exo_list = []
    df_testing_list = []
    df_training_list = []
    RCs_list = []

    for sample_name in df_ann["Sample Name"].tolist():
        df1, df2, df3, df4, dict5 = predict_per_sample(df_cryp, df_exo, df_cs,
                                                       sample_name)
        df_cryp_list.append(df1)
        df_exo_list.append(df2)
        df_testing_list.append(df3)
        df_training_list.append(df4)
        RCs_list.append(dict5)

    df_cryp = pd.concat(df_cryp_list, ignore_index=True)
    df_exo = pd.concat(df_exo_list, ignore_index=True)
    df_testing = pd.concat(df_testing_list, ignore_index=True)
    df_training = pd.concat(df_training_list, ignore_index=True)

    RCs = {}
    for RC in RCs_list:
        for k, v in RC.items():
            RCs[k] = v

    exo_score_thresholds = {}
    cryp_score_thresholds = {}
    with PdfPages(output_pdf) as pdf:
        for sample in df_ann["Sample Name"].tolist():
            df_testing_sample = df_testing[df_testing["Sample Name"] == sample]
            df_exo_sample = df_exo[df_exo["Sample Name"] == sample]
            df_cryp_sample = df_cryp[df_cryp["Sample Name"] == sample]
            if df_testing_sample.shape[0] < 10 or df_exo_sample.shape[
                    0] < 10 or df_cryp_sample.shape[0] < 10:
                continue
            fig, ax, exo_score_threshold, cryp_score_threshold = plot(
                df_testing_sample, df_exo_sample, df_cryp_sample)
            if exo_score_threshold >= 0:
                exo_score_thresholds[sample] = exo_score_threshold
            if cryp_score_threshold >= 0:
                cryp_score_thresholds[sample] = cryp_score_threshold
            fig.suptitle(sample, fontsize=8)
            pdf.savefig(fig)
            plt.close(fig)

    df_exo["RT threshold"] = df_exo["Sample Name"].map(
        lambda x: exo_score_thresholds.get(x, np.nan))
    df_cryp["RT threshold"] = df_cryp["Sample Name"].map(
        lambda x: cryp_score_thresholds.get(x, np.nan))

    pkl.dump(RCs, open(RCs_output, "wb"))
    cols = ["spectrum", "peptide", "RT exp", "Sample Name", "pride id"]
    df_training[cols].to_csv(df_training_output, sep="\t", header=True)
    cols = [
        "spectrum", "peptide", "RT exp", "RT pred", "Sample Name", "pride id"
    ]
    df_testing[cols].to_csv(df_testing_output, sep="\t", header=True)
    cols = [
        "spectrum", "peptide", "RT exp", "RT pred", "predicted score",
        "RT threshold", "Sample Name", "pride id"
    ]
    cryp_r2 = r2_score(df_cryp['RT exp'], df_cryp['RT pred'])
    print(f"cryptic peptides r2 score: {cryp_r2}")
    df_cryp[cols].to_csv(df_cryp_output, sep="\t", header=True)
    df_exo[cols].to_csv(df_exo_output, sep="\t", header=True)


if __name__ == "__main__":
    main()
