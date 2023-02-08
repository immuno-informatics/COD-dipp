import math

import numpy as np
import pandas as pd

import pybedtools
from Bio import SeqIO
from Bio.Seq import Seq
from textwrap import wrap
import argparse
import subprocess
import os


def run_command(command):
    p = subprocess.Popen(command,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    return iter(p.stdout.readline, b'')


def genomic_seq(genome_dict, chrom, start, end):
    return genome_dict[chrom][start - 1:end]


def translate_seq(string, strand):
    if strand == "+":
        return str(Seq(string).translate(to_stop=False))
    else:
        return str(Seq(string).reverse_complement().translate(to_stop=False))


def translate_genomic_seq(genome_dict, chrom, start, end, strand):
    genomic_sequence = genomic_seq(genome_dict, chrom, start, end)

    if strand == "+":
        while len(genomic_sequence) % 3 != 0:
            genomic_sequence = genomic_sequence[1:]
    else:
        while len(genomic_sequence) % 3 != 0:
            genomic_sequence = genomic_sequence[:-1]
    return translate_seq(genomic_sequence, strand)


def inframe_genomic_seq(genome_dict, chrom, start, end, strand):
    genomic_sequence = genomic_seq(genome_dict, chrom, start, end)
    while len(genomic_sequence) % 3 != 0:
        if strand == "+":
            genomic_sequence = genomic_sequence[1:]
        else:
            genomic_sequence = genomic_sequence[:-1]
    if strand == "-":
        genomic_sequence = str(Seq(genomic_sequence).reverse_complement())
    return genomic_sequence


# kozak pwm from "An analysis of S'-noncoding sequences from 699 vertebrate messenger RNAs"
def titer_sequence(genome_dict, chrom, position, strand, len_pep):
    start = position - 201
    end = position + 200
    local_position = 201 if strand == "+" else 201 - len_pep * 3
    local_pos_index = int((local_position / 3) - 1)
    sequence = genomic_seq(genome_dict, chrom, start, end)
    if strand == "-":
        sequence = str(Seq(sequence).reverse_complement())
    tis_codon_list = ["ATG", "CTG", "ACG", "TTG", "GTG"]
    stop_codon_list = ["TAG", "TAA", "TGA"]
    assert len(sequence) % 3 == 0, "sequence length should be dividable by 3"
    sequence = wrap(sequence, 3)
    tis_index = []
    for i in range(local_pos_index, -1, -1):
        codon = sequence[i]
        if codon in stop_codon_list:
            break
        if codon in tis_codon_list:
            tis_index.append(abs(i))
    sequences = []
    sequence = "".join(sequence)
    if tis_index:
        for index in tis_index:
            pos = index * 3
            global_pos = start + pos - 3 if strand == "+" else end - pos + 1
            loc_start = pos - 100 if pos - 100 >= 0 else 0
            loc_end = pos + 3 + 100 if pos + 3 + 100 <= len(sequence) else len(
                sequence)
            padding_start = "$" * abs(100 - (pos - loc_start))
            padding_end = "$" * abs(100 - (loc_end - pos - 3))
            sequences.append([
                padding_start + sequence[loc_start:loc_end] + padding_end,
                global_pos
            ])
    return sequences


def inframe_fun(row):
    if row.strand_coords == "+":
        cond_inframe = (row.start_coords -
                        row["inframe positions start"]) % 3 == 0
    else:
        cond_inframe = (row["inframe positions start"] -
                        row.end_coords) % 3 == 0
    return cond_inframe


def check_methionine(x, aminoacid):
    peptide = x.database_peptide_coords
    intronc_seq = x["intronic seq"]
    try:
        start = intronc_seq.index(peptide)
    except ValueError:
        return np.nan
    seq = intronc_seq[:start][::-1]
    m = 0
    s = 0
    for i in seq:
        if i == "*":
            s += 1
            break
        if i == aminoacid:
            m += 1
    return m > 0


def args():
    parser = argparse.ArgumentParser(
        description=
        "prepares denovo peptides for a TIS analysis using TITER (doi: 10.1093/bioinformatics/btx247)"
    )
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-input_denovo',
        type=str,
        help='Path to tsv containing denovo peptides coordinates annotation',
        required=True)

    requiredNamed.add_argument('-input_genome',
                               type=str,
                               help='Path to hg38 genome fasta',
                               required=True)

    requiredNamed.add_argument(
        '-input_gene_frames',
        type=str,
        help=
        'Path to tsv with gene exonic features start and end inframe coordinates',
        required=True)

    requiredNamed.add_argument(
        '-output_dir',
        type=str,
        help=
        'Path to output dir (outputs: titer_input.txt, df_titer_pos.tsv, df_titer_neg.tsv, 3FT_coords_annotation_3m_inframe.tsv)',
        required=True)

    return parser


if __name__ == "__main__":
    arguments = args().parse_args()
    file_3ft = arguments.input_denovo
    genome_file = arguments.input_genome
    file_frames = arguments.input_gene_frames
    output_dir = arguments.output_dir

    genome_dict = {}
    for rec in SeqIO.parse(genome_file, "fasta"):
        genome_dict[rec.name] = str(rec.seq)

    df_frames = pd.read_csv(file_frames, sep="\t")
    df_3m = pd.read_csv(file_3ft, sep="\t")

    df_frames.chrom = df_frames.chrom.str.replace("chr", "")
    df_frames["strand"] = df_frames["Intron end"].map(
        lambda x: "+" if not math.isnan(x) else "-")
    cols_pos = [
        "chrom", "start", "Intron end", "strand", "Intron start", "end",
        "inframe positions start", "inframe positions end", "name"
    ]
    cols_neg = [
        "chrom", "Intron start", "end", "strand", "Intron end", "start",
        "inframe positions start", "inframe positions end", "name"
    ]

    df_frames["chrom"] = df_frames["chrom"].str.replace("X", "23")
    df_frames["chrom"] = df_frames["chrom"].str.replace("Y", "24")
    df_frames = df_frames[(~df_frames["inframe positions start"].isnull())
                          & (~df_frames["inframe positions end"].isnull())]

    df_frames_neg = df_frames[df_frames["strand"] ==
                              "-"][cols_neg].reset_index(drop=True).copy(
                                  deep=True)
    df_frames_neg["Intron start"] = df_frames_neg["Intron start"].map(int)

    df_frames_pos = df_frames[df_frames["strand"] ==
                              "+"][cols_pos].reset_index(drop=True).copy(
                                  deep=True)
    df_frames_pos["Intron end"] = df_frames_pos["Intron end"].map(int)

    dup_cols = [
        "chr_coords", "start_coords", "end_coords", "strand_coords",
        "feature_ann"
    ]
    df_list = [
        df_3m[df_3m.feature_ann != "Intron"],
        df_3m[df_3m.feature_ann == "Intron"]
    ]
    df_introns = pd.concat(df_list, axis=0).drop_duplicates(dup_cols)
    df_introns = df_introns[df_introns.feature_ann == "Intron"]

    cols_introns = [
        "chr_coords", "start_coords", "end_coords", "strand_coords",
        "denovo_peptide_coords", "database_peptide_coords", "mismatch_coords",
        "gene_coords", "name_ann"
    ]
    df_introns = df_introns[cols_introns]

    df_introns["chr_coords"] = df_introns["chr_coords"].map(str).str.replace(
        "X", "23")
    df_introns["chr_coords"] = df_introns["chr_coords"].map(str).str.replace(
        "Y", "24")

    bed_ann_pos = pybedtools.BedTool.from_dataframe(df_frames_pos)
    bed_ann_neg = pybedtools.BedTool.from_dataframe(df_frames_neg)
    bed_coords = pybedtools.BedTool.from_dataframe(df_introns)
    bed_inter_pos = bed_coords.intersect(bed_ann_pos, wo=True)
    bed_inter_neg = bed_coords.intersect(bed_ann_neg, wo=True)
    # bed_inter = bed_coords.intersect(bed_ann, loj=True)

    cols = []
    cols = df_introns.columns.tolist()
    cols.extend(df_frames_pos.columns.tolist())
    cols.append("overlap")

    df_inter_pos = bed_inter_pos.to_dataframe(disable_auto_names=True)
    df_inter_pos = df_inter_pos.T.reset_index().T.drop_duplicates()
    df_inter_pos.columns = cols
    df_inter_pos = df_inter_pos[df_inter_pos['name'] ==
                                df_inter_pos['name_ann']]
    df_inter_pos = df_inter_pos[df_inter_pos.strand_coords ==
                                df_inter_pos.strand]
    df_inter_pos["chr_coords"] = df_inter_pos["chr_coords"].replace(23, "X")
    df_inter_pos["chr_coords"] = df_inter_pos["chr_coords"].replace(24, "Y")

    cols = []
    cols = df_introns.columns.tolist()
    cols.extend(df_frames_neg.columns.tolist())
    cols.append("overlap")

    df_inter_neg = bed_inter_neg.to_dataframe(disable_auto_names=True)
    df_inter_neg = df_inter_neg.T.reset_index().T.drop_duplicates()
    df_inter_neg.columns = cols
    df_inter_neg = df_inter_neg[df_inter_neg['name'] ==
                                df_inter_neg['name_ann']]
    df_inter_neg = df_inter_neg[df_inter_neg.strand_coords ==
                                df_inter_neg.strand]
    df_inter_neg["chr_coords"] = df_inter_neg["chr_coords"].replace(23, "X")
    df_inter_neg["chr_coords"] = df_inter_neg["chr_coords"].replace(24, "Y")

    df_inter_pos["inframe intron"] = df_inter_pos.apply(inframe_fun, axis=1)
    df_inter_neg["inframe intron"] = df_inter_neg.apply(inframe_fun, axis=1)

    # positive strand methionine check
    # check CUG (L) for codon start
    df_inter_pos["intronic seq"] = df_inter_pos.apply(
        lambda row: translate_genomic_seq(genome_dict, str(row.chr_coords), row
                                          .start, row["end_coords"], "+"),
        axis=1)

    df_inter_pos["genomic seq"] = df_inter_pos.apply(
        lambda row: inframe_genomic_seq(genome_dict, str(row.chr_coords), row.
                                        start, row["end_coords"], "+"),
        axis=1)

    df_inter_pos["inframe intron Met"] = df_inter_pos.apply(
        lambda x: check_methionine(x, "M"), axis=1)

    # negative strand methionine check
    df_inter_neg["intronic seq"] = df_inter_neg.apply(
        lambda row: translate_genomic_seq(genome_dict, str(row.chr_coords),
                                          row["start_coords"], row.end, "-"),
        axis=1)

    df_inter_neg["genomic seq"] = df_inter_neg.apply(
        lambda row: inframe_genomic_seq(genome_dict, str(row.chr_coords), row[
            "start_coords"], row.end, "-"),
        axis=1)

    df_inter_neg["inframe intron Met"] = df_inter_neg.apply(
        lambda x: check_methionine(x, "M"), axis=1)

    met_pep_pos = df_inter_pos.groupby("database_peptide_coords").apply(
        lambda x: any(x["inframe intron Met"]))
    met_pep_neg = df_inter_neg.groupby("database_peptide_coords").apply(
        lambda x: any(x["inframe intron Met"]))

    df = pd.concat([df_inter_pos, df_inter_neg], ignore_index=True)

    agg_params = (("inframe intron", lambda x: any(x["inframe intron"])),
                  ("inframe intron Met",
                   lambda x: any(x["inframe intron Met"])))

    temp1 = df.groupby(["gene_coords",
                        "denovo_peptide_coords"])["inframe intron"].agg([
                            ("exon inframe", lambda x: any(x))
                        ]).reset_index()
    temp2 = df.groupby(["gene_coords",
                        "denovo_peptide_coords"])["inframe intron Met"].agg([
                            ("methionine inframe", lambda x: any(x))
                        ]).reset_index()
    df_final = temp1.merge(temp2, how="outer")

    df_inter_pos["titer seq"] = df_inter_pos.apply(
        lambda row: titer_sequence(genome_dict, str(row.chr_coords), row[
            "start_coords"], "+", len(row.denovo_peptide_coords)),
        axis=1)
    df_titer_pos = df_inter_pos[
        df_inter_pos["titer seq"].map(len) > 0].explode("titer seq")[[
            "chr_coords", "start_coords", "end_coords", "strand_coords",
            "denovo_peptide_coords", "titer seq"
        ]]
    df_titer_pos["tis pos"] = df_titer_pos["titer seq"].map(lambda x: x[1])
    df_titer_pos["titer seq"] = df_titer_pos["titer seq"].map(lambda x: x[0])

    df_inter_neg["titer seq"] = df_inter_neg.apply(
        lambda row: titer_sequence(genome_dict, str(row.chr_coords), row[
            "start_coords"], "-", len(row.denovo_peptide_coords)),
        axis=1)
    df_titer_neg = df_inter_neg[
        df_inter_neg["titer seq"].map(len) > 0].explode("titer seq")[[
            "chr_coords", "start_coords", "end_coords", "strand_coords",
            "denovo_peptide_coords", "titer seq"
        ]]
    df_titer_neg["tis pos"] = df_titer_neg["titer seq"].map(lambda x: x[1])
    df_titer_neg["titer seq"] = df_titer_neg["titer seq"].map(lambda x: x[0])

    titer_seqs = pd.concat([
        df_titer_neg["titer seq"].drop_duplicates(),
        df_titer_pos["titer seq"].drop_duplicates()
    ],
                           ignore_index=True).tolist()

    with open(os.path.join(output_dir, "titer_input.txt"), "w") as out:
        out.write("\n".join(titer_seqs))

    df.to_csv(os.path.join(output_dir, "3ft_coords_annotation_3m_inframe.tsv"),
              sep="\t",
              header=True,
              index=False)
    df_titer_pos.to_csv(os.path.join(output_dir, "df_titer_pos.tsv"),
                        sep="\t",
                        header=True,
                        index=False)
    df_titer_neg.to_csv(os.path.join(output_dir, "df_titer_neg.tsv"),
                        sep="\t",
                        header=True,
                        index=False)
