import argparse

import pandas as pd


def args():
    parser = argparse.ArgumentParser(description='Filters Blat output')
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument('-input_tsv',
                               type=str,
                               help='Path to tsv of denovo BLAT hits',
                               required=True)

    requiredNamed.add_argument(
        '-nunique_filter',
        type=int,
        help='The least number of unique amino acids per sequence (inclusive)',
        required=True)

    requiredNamed.add_argument('-output_tsv',
                               type=str,
                               help='Path to filtered tsv output',
                               required=True)

    return parser


if __name__ == "__main__":
    arguments = args().parse_args()
    input_tsv = arguments.input_tsv
    output_tsv = arguments.output_tsv
    nunique_filter = arguments.nunique_filter

    # loading dataframes
    df = pd.read_csv(input_tsv, sep="\t", header=0)
    cond = df.denovo_seq_nomod.map(lambda x: len(set(x))) >= nunique_filter
    df[cond].to_csv(output_tsv, sep="\t", header=True, index=False)
