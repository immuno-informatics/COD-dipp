import argparse
import math
import os

import numpy as np
import pandas as pd
from cruzdb import Genome


def args():
    parser = argparse.ArgumentParser(
        description='Generates UCSC KnownGene Table, '
        'convert to 1 based coords and extract features (Exon, Intron, UTR)')
    requiredNamed = parser.add_argument_group('required named arguments')

    requiredNamed.add_argument(
        '-output_dir', type=str, help='path to output dir', required=True)

    return parser


arguments = args().parse_args()
outdir = arguments.output_dir

tsv_file1 = os.path.join(outdir, "UCSC_knownGene_hg38.tsv")
tsv_file2 = os.path.join(outdir, "UCSC_knownGene_hg38_features.tsv")

g = Genome(db="hg38")
df = g.dataframe('knownGene')
df["ensembl_transcript_id"] = df.name.str.rsplit(".", 1, expand=True)[0]
df.to_csv(tsv_file1, sep="\t", header=True, index=False)

# split the dataframe to get an exon per line
sep = ","
tsv_sep = "\t"

lines = open(tsv_file1).readlines()
header = lines[0]

outheader = ["name", "chrom", "strand", "start", "end", "feature"]

cols = {}
for i, h in enumerate(header.strip().split(tsv_sep)):
    cols[h] = i

with open(tsv_file2, "w") as out:
    out.write(tsv_sep.join(outheader) + "\n")
    for line in lines[1:]:
        utr5End = ""
        utr5Start = ""
        utr3Start = ""
        utr3End = ""
        l = line.strip().split(tsv_sep)
        # convert to 1 based coord
        l[cols["txStart"]] = str(int(l[cols["txStart"]]) + 1)
        l[cols["cdsStart"]] = str(int(l[cols["cdsStart"]]) + 1)
        l[cols["exonStarts"]] = ",".join([
            str(int(x) + 1)
            for x in l[cols["exonStarts"]].strip(sep).split(sep)
        ])
        l[cols["name"]] = l[cols["name"]].rsplit(".")[0]
        sl = [el.strip(sep).split(sep) for el in l]
        maxf = max([len(el) for el in sl])
        for i in range(len(sl)):
            while len(sl[i]) < maxf:
                sl[i].append(sl[i][-1])
        output_lines = []
        for i in range(maxf):
            output_lines.append([el[i] for el in sl])
        txStart = output_lines[0][cols["txStart"]]
        txEnd = output_lines[0][cols["txEnd"]]
        cdsStart = output_lines[0][cols["cdsStart"]]
        cdsEnd = output_lines[0][cols["cdsEnd"]]
        strand = output_lines[0][cols["strand"]]
        if int(cdsStart) == int(cdsEnd) + 1:
            continue
        if strand == "+" and txStart != cdsStart:
            utr5Start = txStart
            utr5End = str(int(cdsStart) - 1)
            utr3Start = str(int(cdsEnd) + 1)
            utr3End = txEnd
        elif strand == "-" and cdsEnd != txEnd:
            utr5End = txEnd
            utr5Start = str(int(cdsEnd) + 1)
            utr3Start = txStart
            utr3End = str(int(cdsStart) - 1)
        else:
            utr5End = ""
            utr5Start = ""
            utr3Start = ""
            utr3End = ""
        if utr5Start:
            new_line = [l[cols[el]] for el in ["name", "chrom", "strand"]]
            new_line.extend([utr5Start, utr5End, "5UTR"])
            out.write(tsv_sep.join(new_line) + "\n")
        if utr3Start:
            new_line = [l[cols[el]] for el in ["name", "chrom", "strand"]]
            new_line.extend([utr3Start, utr3End, "3UTR"])
            out.write(tsv_sep.join(new_line) + "\n")
        nexons = len(output_lines)
        for i, line in enumerate(output_lines):
            if line[cols["strand"]] == "+":
                if utr5End != "":
                    if int(line[cols["exonEnds"]]) <= int(utr5End):
                        continue
                    else:
                        if int(line[cols["exonStarts"]]) < int(utr5End):
                            line[cols["exonStarts"]] = str(int(utr5End) + 1)
                if utr3Start != "":
                    if int(line[cols["exonStarts"]]) >= int(utr3Start):
                        continue
                    else:
                        if int(line[cols["exonEnds"]]) > int(utr3Start):
                            line[cols["exonEnds"]] = str(int(utr3Start) - 1)
            else:
                if utr5Start != "":
                    if int(line[cols["exonStarts"]]) >= int(utr5Start):
                        continue
                    else:
                        if int(line[cols["exonEnds"]]) > int(utr5Start):
                            line[cols["exonEnds"]] = str(int(utr5Start) - 1)
                if utr3End != "":
                    if int(line[cols["exonEnds"]]) <= int(utr3End):
                        continue
                    else:
                        if int(line[cols["exonStarts"]]) < int(utr3End):
                            line[cols["exonStarts"]] = str(int(utr3End) + 1)
            output_line = [
                line[cols[el]] for el in
                ["name", "chrom", "strand", "exonStarts", "exonEnds"]
            ]
            output_line.append("Exon")
            out.write(tsv_sep.join(output_line) + "\n")
            if i < nexons - 1 and nexons > 1:
                intronStart = str(int(line[cols["exonEnds"]]) + 1)
                intronEnd = str(
                    int(output_lines[i + 1][cols["exonStarts"]]) - 1)
                new_line = [l[cols[el]] for el in ["name", "chrom", "strand"]]
                new_line.extend([intronStart, intronEnd, "Intron"])
                out.write(tsv_sep.join(new_line) + "\n")
