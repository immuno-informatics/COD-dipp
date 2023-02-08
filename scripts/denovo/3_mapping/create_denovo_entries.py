import argparse
from collections import defaultdict

import pandas as pd


def args():
    parser = argparse.ArgumentParser(
        description=
        "adds proteins id, start and end to devo peptide using blat output")
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-input_denovo',
        type=str,
        help='Path to denovo peptides tsv file (deepnovo output)',
        required=True)

    requiredNamed.add_argument(
        '-input_blat',
        type=str,
        help='Path to denovo exonic peptide after blat alignment pslx file',
        required=True)

    requiredNamed.add_argument(
        '-output', type=str, help='Path to output tsv', required=True)

    return parser


global C
C = {
    "matches": 0,
    "misMatches": 1,
    "repMatches": 2,
    "nCount": 3,
    "qNumInsert": 4,
    "qBaseInsert": 5,
    "tNumInsert": 6,
    "tBaseInsert": 7,
    "strand": 8,
    "qName": 9,
    "qSize": 10,
    "qStart": 11,
    "qEnd": 12,
    "tName": 13,
    "tSize": 14,
    "tStart": 15,
    "tEnd": 16,
    "blockCount": 17,
    "blockSizes": 18,
    "qStarts": 19,
    "tStarts": 20,
    "qseq": 21,
    "tseq": 22
}


def parse_pslx(file_path):
    fdict = {}
    for key in ["denovo_seq_nomod", "protein", "start", "end"]:
        fdict[key] = []
    line_d = defaultdict(lambda: defaultdict(list))
    with open(file_path, 'r') as fh:
        for line in fh:
            line = line.strip().split("\t")
            seq = line[C["qseq"]].strip(",")
            line_d[seq]["protein"].append(line[C["tName"]])
            line_d[seq]["start"].append(line[C["tStart"]])
            line_d[seq]["end"].append(line[C["tEnd"]])
    for k1 in line_d.keys():
        fdict["denovo_seq_nomod"].append(k1)
        for k2 in ["protein", "start", "end"]:
            value = ";".join(line_d[k1][k2])
            fdict[k2].append(value)
    return fdict


arguments = args().parse_args()
input_denovo = arguments.input_denovo
input_pslx = arguments.input_blat
output = arguments.output

df_dn = pd.read_csv(input_denovo, sep="\t")
df_pslx = pd.DataFrame(parse_pslx(input_pslx))
df_dn = df_dn.merge(df_pslx, how="inner")
df_dn.to_csv(output, sep="\t", header=True, index=False)
