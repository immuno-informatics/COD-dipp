import argparse
import pickle
import sys

from Bio import SeqIO
from mpi4py import MPI


def args():
    parser = argparse.ArgumentParser(description='Filters Blat output')
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument(
        '-input_fasta',
        nargs=1,
        type=str,
        help='Path to input fasta',
        required=True)

    requiredNamed.add_argument(
        '-input_exon_pslx',
        type=str,
        nargs="*",
        help=
        'Path to first iteration pslx against humandb for introns filtering',
        required=True)

    requiredNamed.add_argument(
        '-exonstats_output',
        type=str,
        help='path to exon stat pickle file output',
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
    input_fasta = arguments.input_fasta[0]
    input_exons_pslx = arguments.input_exon_pslx
    exon_stat_output = arguments.exonstats_output

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    if rank == 0:
        print(f"Number of processes: {size}")
        db = {}
        for seq_record in SeqIO.parse(input_fasta, 'fasta'):
            key = seq_record.id.split(';', 1)[0]
            db[key] = len(str(seq_record.seq))
    else:
        db = None

    db = comm.bcast(db, root=0)
    if rank == 0:
        file_list = chunker_list(input_exons_pslx, size - 1)
        for i, file_names in enumerate(file_list):
            comm.send(file_names, dest=i + 1, tag=i + 1)

    if rank != 0:
        file_names = comm.recv(source=0, tag=rank)
        exon_stat = {}
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
                        if exon_stat[qname] > missed_aa:
                            exon_stat[qname] = missed_aa
                    except KeyError:
                        exon_stat[qname] = missed_aa
        comm.send(exon_stat, dest=0, tag=rank)

    if rank == 0:
        output_data = []
        for i in range(1, size):
            output_data.append(comm.recv(source=i, tag=i))
        exon_stats = {}
        for d in output_data:
            for k, v in d.items():
                try:
                    if exon_stats[k] > v:
                        exon_stats[k] = v
                except KeyError:
                    exon_stats[k] = v
        pickle.dump(exon_stats, open(exon_stat_output, "wb"))


if __name__ == "__main__":
    main()
