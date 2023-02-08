import argparse

import pandas as pd
import pybedtools
from Bio import SeqIO


def args():
    parser = argparse.ArgumentParser(description="Maps denovo's blat "
                                     "filtered output to the genome")
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-input_tsv',
        type=str,
        help='Path to TSV of 3FT filtered versus human database',
        required=True)

    requiredNamed.add_argument(
        '-input_annotation',
        type=str,
        help='path to genomic annotation (UCSC_knownGene_hg38_features.tsv)',
        required=True)

    requiredNamed.add_argument(
        '-input_fasta',
        type=str,
        help='Path to input fasta (3FTgenes_coding.fasta)',
        required=True)

    requiredNamed.add_argument(
        '-output_coordinates',
        type=str,
        help='Path to the output TSV of denovo peptides coordinates',
        required=True)
    requiredNamed.add_argument(
        '-output_coordinates_annotation',
        type=str,
        help='Path to the output TSV of denovo peptides coordinates '
        'annotation (Exon, Intron, UTR)',
        required=True)
    return parser


if __name__ == "__main__":
    arguments = args().parse_args()
    intron_db = arguments.input_fasta
    intron_tsv = arguments.input_tsv
    ann_file = arguments.input_annotation
    output_coords_tsv = arguments.output_coordinates
    output_coords_ann_tsv = arguments.output_coordinates_annotation

    db = SeqIO.parse(intron_db, format="fasta")
    df = pd.read_csv(intron_tsv, sep="\t", header=0)

    # header: >seq1:chr1:12008|12058-12178|12228:+:ENSG00000223972:F1
    print("retrieving coordinates from database")
    id = []
    chr = []
    start = []
    end = []
    strand = []
    gene = []
    prot = []
    frame = []
    for seq_obj in db:
        seq = seq_obj.seq
        header = seq_obj.name.split(":")
        id.append(header[0])
        chr.append(header[1])
        start.append(header[2].split("|")[0])
        end.append(header[2].split("|")[-1])
        strand.append(header[3])
        gene.append(header[4])
        prot.append(str(seq))
        frame.append(int(header[-1][1]) - 1)

    df_db = pd.DataFrame({
        'id': id,
        'chr': chr,
        'start': start,
        'end': end,
        'strand': strand,
        'gene': gene,
        'frame': frame,
        'prot': prot
    })

    df_db.start = df_db.start.astype(int)
    df_db.end = df_db.end.astype(int)

    print("retrieving coordinates for peptides")
    global df_dict
    df_dict = {}
    df_dict["denovo_peptide"] = []
    df_dict["database_peptide"] = []
    df_dict["chr"] = []
    df_dict["start"] = []
    df_dict["end"] = []
    df_dict["strand"] = []
    df_dict["gene"] = []
    df_dict["mismatch"] = []
    df_dict["frame"] = []

    def extract_coords(df_loc, df_db_loc):
        denovo_peptide = df_loc.qseq
        peptide = df_loc.tseq
        mismatch = df_loc.misMatches
        seqid = df_loc.tName.split(":", 1)[0]
        sdf = df_db_loc[df_db.id == seqid].iloc[0]
        strand = sdf["strand"]
        try:
            start = sdf.prot.index(peptide) * 3 + sdf["frame"]
        except ValueError:
            return
        df_dict["database_peptide"].append(peptide)
        df_dict["denovo_peptide"].append(denovo_peptide)
        if strand == "+":
            df_dict["start"].append(start + sdf["start"])
            df_dict["end"].append(df_dict["start"][-1] + len(peptide) * 3 - 1)
        else:
            df_dict["end"].append(sdf["end"] - start)
            df_dict["start"].append(df_dict["end"][-1] - len(peptide) * 3 + 1)
        df_dict["chr"].append(sdf["chr"])
        df_dict["gene"].append(sdf["gene"])
        df_dict["strand"].append(strand)
        df_dict["mismatch"].append(mismatch)
        df_dict["frame"].append(sdf["frame"])

    df.apply(lambda x: extract_coords(x, df_db), axis=1)

    df_coords = pd.DataFrame(df_dict).drop_duplicates()
    df_coords.to_csv(output_coords_tsv, sep="\t", header=True, index=False)

    print("retrieving peptidic features")
    df_coords = pd.read_csv(output_coords_tsv, sep="\t")
    df_ann = pd.read_csv(ann_file, sep="\t")

    df_coords["score"] = 1000
    df_ann["score"] = 1000

    df_ann["id"] = list(range(df_ann.shape[0]))
    df_coords["id"] = list(range(df_coords.shape[0]))

    df_coords = df_coords[[
        "chr", "start", "end", "id", "score", "strand", "denovo_peptide",
        "database_peptide", "mismatch", "frame", "gene"
    ]]
    df_ann = df_ann[[
        "chrom", "start", "end", "id", "score", "strand", "feature", "name"
    ]]

    df_ann = df_ann[abs(df_ann["start"] - df_ann["end"]) >= 12]

    df_coords["chr"] = df_coords["chr"].str.replace("^chr", "")
    df_ann["chrom"] = df_ann["chrom"].str.replace("^chr", "")

    df_coords["start"] = df_coords["start"].astype(int)
    df_coords["end"] = df_coords["end"].astype(int)

    df_ann["start"] = df_ann["start"].astype(int)
    df_ann["end"] = df_ann["end"].astype(int)

    bed_ann = pybedtools.BedTool.from_dataframe(df_ann)
    bed_coords = pybedtools.BedTool.from_dataframe(df_coords)

    # loj: Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B.
    # s: match strands
    # f:1 cover the full length of A (peptide)

    bed_inter = bed_coords.intersect(bed_ann, loj=True, s=True, f=1)

    cols = df_coords.columns.map(lambda x: x + "_coords").tolist()
    cols.extend(df_ann.columns.map(lambda x: x + "_ann").tolist())

    df_inter = bed_inter.to_dataframe(names=cols)

    df_inter.to_csv(output_coords_ann_tsv, sep="\t", header=True, index=False)
