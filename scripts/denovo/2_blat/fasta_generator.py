import argparse

import pandas as pd


def args():
    parser = argparse.ArgumentParser(
        description=
        'Takes Denovo and closed search outputs and generates fastas')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '-input_closed_search',
        type=str,
        help='Path to tsv file (MS-GF+) peptides TSV',
        required=True)

    requiredNamed.add_argument(
        '-input_denovo',
        type=str,
        help='Path to tsv file (DeepNovo) ALC 90 results',
        required=True)

    requiredNamed.add_argument(
        '-nsplits',
        type=int,
        help='number of splits for the denovo fasta output',
        required=True)

    requiredNamed.add_argument(
        '-output_closed_search',
        type=str,
        help='path to output tsv',
        required=True)

    requiredNamed.add_argument(
        '-output_denovo', type=str, help='path to output tsv', required=True)
    return parser


def slice_per(source, step):
    return [source[i::step] for i in range(step)]


if __name__ == '__main__':
    arguments = args().parse_args()
    file_cs = arguments.input_closed_search
    file_dn = arguments.input_denovo
    output_file_cs = arguments.output_closed_search
    output_file_dn = arguments.output_denovo
    nsplits = arguments.nsplits

    df = pd.read_csv(
        file_dn,
        sep='\t',
        header=0,
        usecols=['feature_id', 'denovo_seq_nomod'])
    df['feature_id'] = df['feature_id'].str.split(":", expand=True)[0]

    global denovo_output_list
    denovo_output_list = []

    def output_fun_dn(df):
        header = ';'.join(df['feature_id'].tolist())
        denovo_output_list.append(f">{header}\n{df.name}")

    df.groupby('denovo_seq_nomod').apply(output_fun_dn)

    if nsplits == 1:
        with open(output_file_dn, 'w') as out:
            out.write("\n".join(denovo_output_list))
    else:
        denovo_output_lists = slice_per(denovo_output_list, nsplits)
        try:
            prefix, suffix = output_file_dn.rsplit('.', 1)
            suffix = '.' + suffix
        except ValueError:
            prefix = output_file_dn
            suffix = ""
        for index, towrite in enumerate(denovo_output_lists):
            split_name = prefix + str(index + 1) + suffix
            with open(split_name, "w") as out:
                out.write("\n".join(towrite))

    df = pd.read_csv(
        file_cs, sep='\t', header=0, usecols=['peptide', 'spectrum'])

    def output_fun_cs(df, fh):
        header = ';'.join(df['spectrum'].tolist())
        fh.write(f">{header}\n{df.name}\n")

    with open(output_file_cs, "w") as fh:
        df.groupby('peptide').apply(lambda x: output_fun_cs(x, fh))