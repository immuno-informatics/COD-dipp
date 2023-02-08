import argparse
import time

import pandas as pd
import pysam
from Bio import SeqIO
import os

def args():
    parser = argparse.ArgumentParser(description="Maps denovo's blat "
                                     "filtered output to the genome")
    requiredNamed = parser.add_argument_group('required named arguments')
    optionalNamed = parser.add_argument_group('optional named arguments')

    requiredNamed.add_argument(
        '-input_coordinates',
        type=str,
        help='Path to denovo peptides coordinates tsv file',
        required=True)

    requiredNamed.add_argument(
        '-input_denovo',
        type=str,
        help='Path to deepnovo psm file',
        required=True)

    requiredNamed.add_argument(
        '-input_genome',
        type=str,
        help='Path to fasta.gz genome',
        required=True)

    requiredNamed.add_argument(
        '-output', type=str, help='path to output sam', required=True)

    optionalNamed.add_argument(
        '-genomic_position_filter',
        type=int,
        default=1,
        help='Filter peptides aligning to more than the provided number',
        required=False)

    return parser


arguments = args().parse_args()
input_file = arguments.input_coordinates
input_denovo = arguments.input_denovo
genome_file = arguments.input_genome
genomic_position_filter = arguments.genomic_position_filter
output = arguments.output

if not output.endswith(".bam"):
    output += ".bam"

print("reading genome")
global genome_dict
genome_dict = {}

for rec in SeqIO.parse(genome_file, "fasta"):
    genome_dict[rec.name] = str(rec.seq)

print("reading genome done")

print("reading data")
df = pd.read_csv(input_file, sep="\t")
denovo = pd.read_csv(input_denovo, sep="\t")

df = df[["denovo_peptide", "chr", "start", "end", "strand", "frame"]]
df.columns = df.columns.map(lambda x: x + "_coords")
df.columns = df.columns.str.replace("denovo_peptide_coords",
                                    "denovo_seq_nomod")

ngenomic = df.groupby("denovo_seq_nomod").apply(
    lambda x: x.drop_duplicates(["chr_coords", "start_coords", "end_coords"]).shape[0]
).to_dict()
df["ngenomic"] = df.denovo_seq_nomod.map(ngenomic)
df = df[df.ngenomic <= genomic_position_filter].drop_duplicates()
df["dna_seq"] = df.apply(lambda x: genome_dict[x.chr_coords.replace("chr", "")][x.start_coords-1:x.end_coords], axis=1)

denovo = denovo.merge(df, how="inner")
denovo["spectrumid"] = denovo["pride id"] + ":" + denovo[
    "sample"] + ":" + denovo["feature_id"]
print("reading data done")

organism = "homo_sapiens"
genome = "genome"

header = {
    "HD": {
        "VN": "1.0.2",
        "SO": "unknown"
    },
    "SQ": [{
        "SN": "chr1",
        "LN": 248956422,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr2",
        "LN": 242193529,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr3",
        "LN": 198295559,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr4",
        "LN": 190214555,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr5",
        "LN": 181538259,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr6",
        "LN": 170805979,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr7",
        "LN": 159345973,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr8",
        "LN": 145138636,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr9",
        "LN": 138394717,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr10",
        "LN": 133797422,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr11",
        "LN": 135086622,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr12",
        "LN": 133275309,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr13",
        "LN": 114364328,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr14",
        "LN": 107043718,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr15",
        "LN": 101991189,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr16",
        "LN": 90338345,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr17",
        "LN": 83257441,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr18",
        "LN": 80373285,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr19",
        "LN": 58617616,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr20",
        "LN": 64444167,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr21",
        "LN": 46709983,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chr22",
        "LN": 50818468,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chrX",
        "LN": 156040895,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chrY",
        "LN": 57227415,
        "AS": genome,
        "SP": organism
    }, {
        "SN": "chrMT",
        "LN": 16569,
        "AS": genome,
        "SP": organism
    }],
    "RG": [{
        "ID": "probam"
    }],
    "CO": ["AS: ENSEMBL", "VN: 94"]
}


def write_alig(row):
    a = pysam.AlignedSegment()
    a.query_name = row.spectrumid
    a.flag = 0x10 if row.strand_coords == "-" else 0x100
    try:
        chrid = int(row.chr_coords.replace('chr', "")) - 1
    except ValueError:
        chrid = 23 if "X" in row.chr_coords else 24
    a.reference_id = chrid
    a.reference_start = row.start_coords - 1 
    a.query_sequence = row.dna_seq
    a.mapping_quality = 255
    a.cigar = [(0, len(row.denovo_seq_nomod) * 3)]
    a.tags = [("NH", row.ngenomic), ("XA", 2), ("XB", "*"),
              ("XC", row.precursor_charge), ("XE", -1), ("XF",
                                                         row.frame_coords),
              ("XG", 0), ("XI", -1), ("XL", -1), ("XM", "*"), ("XN", -1),
              ("XO", "*"), ("XP", row.denovo_seq_nomod), ("XQ", -1), ("XR", "*"),
              ("XS", row.predicted_score_max), ("XU", "*"), ("YA", "*"),
              ("YB", "*"), ("YP", "*")]
    outf.write(a)


print("writing bam file(s)")
global outf

for group, sdf in denovo.groupby(["pride id", "sample"]):
    basename = os.path.basename(output)
    basename = "_".join(list(group)) + "_" + basename
    dirname  = os.path.dirname(output)
    output_local = os.path.join(dirname, basename)
    outf = pysam.AlignmentFile( output_local, "w", header=header)
    denovo.apply(write_alig, axis=1)
    outf.close()
    # print("sorting bam file")
    sorted_output = output_local[:-4] + ".sorted.bam"
    pysam.sort("-o", (sorted_output), (output_local))
    pysam.index(sorted_output)
    # print("writing sam file")
    infile = pysam.AlignmentFile(sorted_output, "rb")
    outfile = pysam.AlignmentFile(sorted_output[:-4] + ".sam", "w", template=infile)
    for s in infile:
        outfile.write(s)
    outfile.close()