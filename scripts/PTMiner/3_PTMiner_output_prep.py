import argparse
import os
import pickle as pkl
from collections import defaultdict
from io import StringIO

import pandas as pd

from Bio import SeqIO


def args():
    parser = argparse.ArgumentParser(
        description=
        "Re-converts PTMiner 1.1.2 output's numerical IDs to Protein Accessions"
    )
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-fasta', type=str, help='Path to input fasta', required=True)

    requiredNamed.add_argument(
        '-db_headers',
        type=str,
        help='Path to the database {"Numerical ID": "Protein Accession"} dict',
        required=True)

    requiredNamed.add_argument(
        '-spectra',
        type=str,
        help=
        'Path to the spectra {"modified spectrum name": "original spectrum name"} dict',
        required=True)

    requiredNamed.add_argument(
        '-filtered',
        type=str,
        help='Path to the "Filtered" results of PTMiner',
        required=True)

    requiredNamed.add_argument(
        '-loc',
        type=str,
        help='Path to the "Localization" results of PTMiner',
        required=True)

    requiredNamed.add_argument(
        '-anno',
        type=str,
        help='Path to the "Annotation" results of PTMiner',
        required=True)

    requiredNamed.add_argument(
        '-output_dir',
        type=str,
        help='Path to output directory location',
        required=True)

    requiredNamed.add_argument(
        '-suffix',
        type=str,
        help='suffix to add for the modifide pepXML input(s) '\
        'before the .txt extension',
        required=True)

    return parser


def read_fasta(fasta_file):
    dict_seq = {}
    for bioseq in SeqIO.parse(fasta_file, format="fasta"):
        dict_seq[bioseq.name] = str(bioseq.seq)
    return dict_seq


def get_prot_pos(x, db):
    pos = x["Position"]
    if pos <= 0:
        pos = 1
    return db[x["Protein Access"].split(" ")[0]].index(x["Sequence"]) + pos


def read_annot(input_file):
    df = pd.read_csv(input_file, sep="\t", header=0)
    header1 = df.columns.tolist()
    header2 = df.iloc[0].dropna().tolist()
    df = df.iloc[1:]

    header2.insert(0, "merger")
    header2[1] = "#"

    ids = df["#"].tolist()
    prev = ""
    for i, id in enumerate(ids):
        if id != "*":
            prev = id
        else:
            ids[i] = prev

    df.insert(0, "merger", ids)

    df1 = df[df["#"] != "*"]

    df2 = df[df["#"] == "*"]
    df2.columns = range(df2.shape[1])
    df2 = df2.loc[:, range(len(header2))]
    df2.columns = header2
    df2.drop("#", axis=1, inplace=True)

    df3 = df1.merge(df2, left_on="merger", right_on="merger", how="left")
    df3.drop("merger", axis=1, inplace=True)
    return df3


if __name__ == "__main__":
    arguments = args().parse_args()
    database_file = arguments.fasta
    file_fil = arguments.filtered
    file_loc = arguments.loc
    file_anno = arguments.anno
    output_dir = arguments.output_dir
    suffix = arguments.suffix
    dict_file = arguments.db_headers
    spec_dict_file = arguments.spectra

    df_fil = pd.read_csv(file_fil, sep="\t")
    df_loc = pd.read_csv(file_loc, sep="\t")
    dict_con = pkl.load(open(dict_file, "rb"))
    db = read_fasta(database_file)
    spectrum_conv = pkl.load(open(spec_dict_file, "rb"))

    try:
        df_ann = read_annot(file_anno)
    except:
        print("reading annotation failed")
        df_ann = None

    df_fil["Protein Access"] = df_fil["Protein Access"].map(dict_con)
    df_loc["Protein Access"] = df_loc["Protein Access"].map(dict_con)
    df_loc["Protein Position"] = df_loc.apply(
        lambda x: get_prot_pos(x, db), axis=1)

    df_loc["Spectrum Name"] = df_loc["Spectrum Name"].map(lambda x: spectrum_conv.get(x, x))
    df_fil["Spectrum Name"] = df_fil["Spectrum Name"].map(lambda x: spectrum_conv.get(x, x))

    df_loc["Dataset Name"] = df_loc["Dataset Name"].map(lambda x: spectrum_conv.get(x, x))
    df_fil["Dataset Name"] = df_fil["Dataset Name"].map(lambda x: spectrum_conv.get(x, x))

    if isinstance(df_ann, pd.DataFrame):
        df_ann["Protein Access"] = df_ann["Protein Access"].map(dict_con)
        df_ann["Protein Position"] = df_ann.apply(
            lambda x: get_prot_pos(x, db), axis=1)
        df_ann["Spectrum Name"] = df_ann["Spectrum Name"].map(lambda x: spectrum_conv.get(x, x))
        df_ann["Dataset Name"] = df_ann["Dataset Name"].map(lambda x: spectrum_conv.get(x, x))
        df_list = [df_fil, df_loc, df_ann]
        file_list = [file_fil, file_loc, file_anno]
    else:
        df_list = [df_fil, df_loc]
        file_list = [file_fil, file_loc]

    for df, file_path in zip(df_list, file_list):
        output = os.path.basename(file_path)
        output = output.rsplit(".", 1)
        output = output[0] + suffix + ".tsv"
        output = os.path.join(output_dir,
                              file_path.rsplit(".")[0] + suffix + ".tsv")
        df.to_csv(output, sep="\t", header=True, index=False)
