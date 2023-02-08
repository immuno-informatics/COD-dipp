import argparse
import sys

from Bio import SeqIO


def args():
    parser = argparse.ArgumentParser(
        description=
        'generate fasta with exlusion or inlucsion of entries identified in pslx'
    )
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument('-filter_type',
                               type=str,
                               choices=["include", "exclude"],
                               required=True)

    requiredNamed.add_argument('-pslx',
                               type=str,
                               help='Path to pslx file',
                               required=True)

    requiredNamed.add_argument('-db_file',
                               type=str,
                               help='path to database file',
                               required=True)

    requiredNamed.add_argument('-output',
                               type=str,
                               help='path to output fasta file',
                               required=True)

    requiredNamed.add_argument(
        '-nsplits',
        type=int,
        help='number of splits for the denovo fasta output',
        required=True)

    return parser


def slice_per(source, step):
    return [source[i::step] for i in range(step)]


if __name__ == '__main__':
    arguments = args().parse_args()
    filter_type = arguments.filter_type
    pslx = arguments.pslx
    db_file = arguments.db_file
    output = arguments.output
    nsplits = arguments.nsplits

    if filter_type == "include":
        cond = [True, False]
    elif filter_type == "exclude":
        cond = [False, True]
    else:
        print("First argument accepted values: include, exclude")
        exit(1)

    qnames = {}
    with open(pslx) as fh:
        for line in fh:
            qnames[line.strip().split("\t")[9].split(";")[0]] = cond[0]

    sequences = []
    with open(db_file) as fasta_file:
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            if qnames.get(seq_record.id.split(';')[0], cond[1]):
                sequences.append(f">{seq_record.id}\n{str(seq_record.seq)}")

    if nsplits == 1:
        with open(output, "w") as out:
            out.write("\n".join(sequences))
    else:
        sequences_splits = slice_per(sequences, nsplits)
        try:
            prefix, suffix = output.rsplit('.', 1)
            suffix = '.' + suffix
        except ValueError:
            prefix = output
            suffix = ""

        for index, towrite in enumerate(sequences_splits):
            if towrite:
                split_name = prefix + str(index + 1) + suffix
                with open(split_name, "w") as out:
                    out.write("\n".join(towrite))
