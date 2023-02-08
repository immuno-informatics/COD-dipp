import argparse
import pickle

from Bio import SeqIO
from mpi4py import MPI


def args():
    parser = argparse.ArgumentParser(description='Filters Blat output')
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-input_pslx',
        type=str,
        nargs="*",
        help='Path to pslx file(s)',
        required=True)

    requiredNamed.add_argument(
        '-exonstats_input',
        type=str,
        help='path to exon stat pickle file input',
        required=True)

    requiredNamed.add_argument(
        '-input_fasta', type=str, help='Path to input fasta', required=True)

    requiredNamed.add_argument(
        '-NmisMatch',
        type=int,
        help='Number of allowed misMatches',
        required=True)

    requiredNamed.add_argument(
        '-NmisMatch_diff',
        type=int,
        default=4,
        help=
        'Filter based on Nmatch(non-exon) - Nmatch(exon) >= NmisMatch_diff',
        required=True)

    requiredNamed.add_argument(
        '-output_pslx',
        type=str,
        help='Path to the filtered pslx output file',
        required=True)

    return parser


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


def chunker_list(seq, size):
    return [seq[i::size] for i in range(size)]


def main():
    arguments = args().parse_args()
    input_pslx_files = arguments.input_pslx
    input_fasta = arguments.input_fasta
    output_pslx = arguments.output_pslx
    NmisMatch = arguments.NmisMatch
    NmisMatch_diff = arguments.NmisMatch_diff
    exon_stat_input = arguments.exonstats_input

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        exon_stat = pickle.load(open(exon_stat_input, "rb"))
        db = {}
        for seq_record in SeqIO.parse(input_fasta, 'fasta'):
            key = seq_record.id.split(';', 1)[0]
            db[key] = len(str(seq_record.seq))
    else:
        db = None
        exon_stat = None

    db = comm.bcast(db, root=0)
    exon_stat = comm.bcast(exon_stat, root=0)

    if rank == 0:
        file_list = chunker_list(input_pslx_files, size - 1)
        for i, file_names in enumerate(file_list):
            comm.send(file_names, dest=i + 1, tag=i + 1)

    if rank != 0:
        file_names = comm.recv(source=0, tag=rank)
        out = []
        for file_name in file_names:
            with open(file_name) as fh:
                for line in fh:
                    line = line.strip().split("\t")
                    blockCount = int(line[C["blockCount"]])
                    if blockCount > 1:
                        continue
                    qname = line[C["qName"]].split(";", 1)[0]
                    qsize = db[qname]
                    matches = int(line[C["matches"]])
                    missed_aa = qsize - matches
                    try:
                        if missed_aa <= NmisMatch and exon_stat[qname] >= NmisMatch_diff:
                            out.append("\t".join(line))
                    except KeyError:
                        if missed_aa <= NmisMatch:
                            out.append("\t".join(line))

        comm.send(out, dest=0, tag=rank)

    if rank == 0:
        with open(output_pslx, "w") as outfile:
            for i in range(1, size):
                output_data = comm.recv(source=i, tag=i)
                outfile.write("\n".join(output_data) + "\n")
                del output_data


if __name__ == "__main__":
    main()
