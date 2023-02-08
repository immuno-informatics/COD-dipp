import argparse
import numpy as np
import pandas as pd
from Bio import SeqIO


def args():
    parser = argparse.ArgumentParser(description='Filters Blat output')
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument('-input_fasta',
                               type=str,
                               help='Path to input fasta',
                               required=True)

    requiredNamed.add_argument(
        '-input_denovoexons',
        type=str,
        help='Path to pslx of denovo versus human database',
        required=True)

    requiredNamed.add_argument(
        '-input_denovo_nonexons',
        type=str,
        help='Path to pslx of denovo versus human alternative frame database',
        required=True)

    requiredNamed.add_argument('-input_denovo90ALC',
                               type=str,
                               help='Path to denovo 90 ALC TSV file',
                               required=True)

    requiredNamed.add_argument('-output_denovoexons',
                               type=str,
                               help='Path to processed denovo exons tsv file',
                               required=True)

    requiredNamed.add_argument(
        '-output_denovo_nonexons',
        type=str,
        help='Path to processed denovo nonexons tsv file',
        required=True)

    return parser


def load_db(input_fasta):
    db = {}
    for seq_record in SeqIO.parse(input_fasta, 'fasta'):
        key = seq_record.id.split(';', 1)[0]
        db[key] = {"header": seq_record.id, "protein": str(seq_record.seq)}
    return db


def retrieve_header(blat_qName, db_dict):
    return db_dict[blat_qName.split(";", 1)[0]]["header"]


def isInDb(seq, db_dict, blat_qName):
    """
        checks if the blat reported qseq exists in the database
        seq: (str) peptide seq
        key: (str) blat reported qName
        db_dict: (dict) output of load_db()
    """
    # blat trims the header if it's longer than a certain threshold.
    key = blat_qName.split(";", 1)[0]
    if seq in db_dict[key]["protein"]:
        return True
    return False


def pos_fun(x):
    return [i for i, (a, b) in enumerate(zip(x.qseq, x.tseq)) if a != b]


def get_pos_scores(x):
    try:
        return ", ".join([str(x[pos]) for pos in x["misMatches pos"]])
    except KeyError:
        return ""


if __name__ == "__main__":
    arguments = args().parse_args()
    input_fasta = arguments.input_fasta
    input_dn_exons = arguments.input_denovoexons
    input_dn_introns = arguments.input_denovo_nonexons
    file_dn = arguments.input_denovo90ALC
    output_dn_exons = arguments.output_denovoexons
    output_dn_introns = arguments.output_denovo_nonexons

    # loading Database
    print("loading Database")
    db = load_db(input_fasta)

    cols = [
        "matches", "misMatches", "repMatches", "nCount", "qNumInsert",
        "qBaseInsert", "tNumInsert", "tBaseInsert", "strand", "qName", "qSize",
        "qStart", "qEnd", "tName", "tSize", "tStart", "tEnd", "blockCount",
        "blockSizes", "qStarts", "tStarts", "qseq", "tseq"
    ]

    # loading dataframes
    print("loading dataframes")
    df_dn_exons = pd.read_csv(input_dn_exons,
                              sep="\t",
                              header=None,
                              names=cols)

    df_dn_introns = pd.read_csv(input_dn_introns,
                                sep="\t",
                                header=None,
                                names=cols)

    df_dn = pd.read_csv(file_dn,
                        sep='\t',
                        header=0,
                        usecols=[
                            "feature_id", "predicted_score",
                            "predicted_score_max", "predicted_position_score",
                            "precursor_charge", "precursor_mz",
                            "denovo_seq_nomod"
                        ])

    print("converting feature id to spectrum id")
    df_dn['qName'] = df_dn['feature_id'].str.split(":", expand=True)[0]

    print("scaling denovo score between 0-100")
    df_dn["predicted_score"] = 100 * np.exp(df_dn["predicted_score"])

    # discarding BLAT hits with mismatches
    print("Keeping hits with 0 mismatch (I and L exception)")

    df_dn_introns = df_dn_introns[df_dn_introns.qseq.str.replace("I", "L") ==
                                  df_dn_introns.tseq.str.replace("I", "L")]
    df_dn_exons = df_dn_exons[df_dn_exons.qseq.str.replace("I", "L") ==
                              df_dn_exons.tseq.str.replace("I", "L")]

    df_dn_introns.tseq = df_dn_introns.tseq.map(lambda x: x.strip(','))
    df_dn_introns.qseq = df_dn_introns.qseq.map(lambda x: x.strip(','))

    df_dn_exons.tseq = df_dn_exons.tseq.map(lambda x: x.strip(','))
    df_dn_exons.qseq = df_dn_exons.qseq.map(lambda x: x.strip(','))

    # Removing sequences that do not exists in the database
    # side effect of having '*' representing a stop codon in the translated
    # sequence. Blat removes these characters and joins the sequences.
    # Hence some peptides might not exist in the database
    cond = df_dn_introns.apply(lambda x: isInDb(x.qseq, db, x.qName), axis=1)
    df_dn_introns = df_dn_introns[cond]

    df_dn_exons["isDecoy"] = df_dn_exons.tName.str.contains("^rev")

    df_dn_exons["qName id"] = df_dn_exons.qName.astype("category").cat.codes
    df_dn_introns["qName id"] = df_dn_introns.qName.astype(
        "category").cat.codes

    print("dropping duplicates for mapping (for speed)")
    df_dn_exons_reduced = df_dn_exons[["qName", "tseq", "qseq",
                                       "qName id"]].drop_duplicates()
    df_dn_introns_reduced = df_dn_introns[[
        "qName", "tseq", "qseq", "qName id"
    ]].drop_duplicates()

    # uncomment if you want to consider BLAT hits with mismatches

    # print("retrieving the mismatch positions")
    # df_dn_introns_reduced["misMatches pos"] = df_dn_introns_reduced.apply(
    #     pos_fun, axis=1)
    # df_dn_exons_reduced["misMatches pos"] = df_dn_exons_reduced.apply(pos_fun,
    #                                                                   axis=1)

    print("regenerating the full qName in case it got trimmed by BLAT")
    df_dn_exons_reduced.qName = df_dn_exons_reduced.apply(
        lambda x: retrieve_header(x.qName, db), axis=1)
    df_dn_introns_reduced.qName = df_dn_introns_reduced.apply(
        lambda x: retrieve_header(x.qName, db), axis=1)

    print("generating tidy format tables")
    df_dn_introns_reduced_tidy = df_dn_introns_reduced.copy(deep=True)
    df_dn_introns_reduced_tidy["qName"] = df_dn_introns_reduced[
        "qName"].str.split(";")
    df_dn_introns_reduced_tidy = df_dn_introns_reduced_tidy.explode("qName")

    df_dn_exons_reduced_tidy = df_dn_exons_reduced.copy(deep=True)
    df_dn_exons_reduced_tidy["qName"] = df_dn_exons_reduced_tidy[
        "qName"].str.split(";")
    df_dn_exons_reduced_tidy = df_dn_exons_reduced_tidy.explode("qName")

    print("adding psm level mass spectrometry information to peptides")
    df_dn_introns_reduced_tidy = df_dn.merge(df_dn_introns_reduced_tidy,
                                             how="right")
    df_dn_exons_reduced_tidy = df_dn.merge(df_dn_exons_reduced_tidy,
                                           how="right")

    # uncomment if you want to consider BLAT hits with mismatches

    # print("scaling denovo score between 0-100")
    # scores_list = [
    #     x.split(",") for x in df_dn["predicted_position_score"].values
    # ]
    # predicted_position_score = 100 * np.exp(
    #     pd.DataFrame(scores_list, dtype="float"))
    # predicted_position_score["qName"] = df_dn["qName"].tolist()

    # print("adding the mismatch positional scores")
    # temp_df = predicted_position_score.merge(
    #     df_dn_introns_reduced_tidy[["qName", "misMatches pos"]])
    # temp_df["misMatches pos scores"] = temp_df.apply(get_pos_scores, axis=1)
    # df_dn_introns_reduced_tidy = df_dn_introns_reduced_tidy.merge(
    #     temp_df[["qName", "misMatches pos scores"]], how="left")

    # temp_df = predicted_position_score.merge(
    #     df_dn_exons_reduced_tidy[["qName", "misMatches pos"]])
    # temp_df["misMatches pos scores"] = temp_df.apply(get_pos_scores, axis=1)
    # df_dn_exons_reduced_tidy = df_dn_exons_reduced_tidy.merge(
    #     temp_df[["qName", "misMatches pos scores"]], how="left")

    # print('converting "misMatches pos" column from type (list) '
    #       'to type (string) with a ", " seperator')
    # df_dn_introns_reduced_tidy["misMatches pos"] = df_dn_introns_reduced_tidy[
    #     "misMatches pos"].map(lambda x: ", ".join([str(y) for y in x]))

    # df_dn_exons_reduced_tidy["misMatches pos"] = df_dn_exons_reduced_tidy[
    #     "misMatches pos"].map(lambda x: ", ".join([str(y) for y in x]))

    print("writing output")
    df_dn_introns_reduced_tidy.merge(df_dn_introns.drop("qName", axis=1),
                                     how="inner").to_csv(output_dn_introns,
                                                         sep="\t",
                                                         header=True,
                                                         index=False)

    # can be too big to fit in memory if considering mismatches
    df_dn_exons_reduced_tidy.merge(df_dn_exons.drop("qName", axis=1),
                                   how="inner").to_csv(output_dn_exons,
                                                       sep="\t",
                                                       header=True,
                                                       index=False)

    # solution: merge chunks
    # creating a empty bucket to save result
    # df1 = df_dn_exons_reduced_tidy
    # df2 = df_dn_exons.drop("qName", axis=1)

    # temp_out = "temp_df.tsv"
    # df2.to_csv(temp_out, header=True, index=False, sep="\t")

    # df_merge = pd.DataFrame(columns=(df1.columns.append(df2.columns)).unique())
    # df_merge.to_csv(output_dn_exons, index=False, sep='\t', header=True)

    # def preprocess(df1, df2, out):
    #     res = pd.merge(df1, df2, how="inner")
    #     res.to_csv(out, mode="a", header=False, index=False, sep='\t')

    # reader = pd.read_csv(temp_out, chunksize=50000, sep="\t", header=0)

    # [preprocess(df1, chunk_df, output_dn_exons) for chunk_df in reader]
    # os.system(f"rm {temp_out}")
