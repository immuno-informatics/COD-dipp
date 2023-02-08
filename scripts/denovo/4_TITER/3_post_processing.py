import argparse
import pandas as pd

def args():
    parser = argparse.ArgumentParser(
        description=
        "post processing for TITER (doi: 10.1093/bioinformatics/btx247)")
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-input_titer',
        type=str,
        help='Path to tsv of TITER output',
        required=True)

    requiredNamed.add_argument(
        '-input_pos',
        type=str,
        help='Path to tsv with positive strand features',
        required=True)

    requiredNamed.add_argument(
        '-input_neg',
        type=str,
        help='Path to tsv with negative strand features',
        required=True)
    requiredNamed.add_argument(
        '-output', type=str, help='Path to output (Titer input)', required=True)

    return parser


if __name__ == "__main__":
    arguments = args().parse_args()
    input_titer = arguments.input_titer
    input_pos = arguments.input_pos
    input_neg = arguments.input_neg
    output = arguments.output

    df_titer_pred = pd.read_csv(input_titer, sep="\t")
    df_titer_pos = pd.read_csv(input_pos, sep="\t")
    df_titer_neg = pd.read_csv(input_neg, sep="\t")
    df_titer_pos = df_titer_pos.merge(df_titer_pred, how="left")
    df_titer_neg = df_titer_neg.merge(df_titer_pred, how="left")
    df_titer_tis = pd.concat([df_titer_pos, df_titer_neg], ignore_index=True)
    df_titer_tis.to_csv(output, sep="\t", index=False, header=True)
