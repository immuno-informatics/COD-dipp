#ecoding:utf-8
"""
    script for converting psm output of philosopher to mzTAB format.
    Author: Georges BEDRAN
"""

import re
import sys

from pyteomics import mass


def get_cols(first_line):
    line = first_line.strip().split('\t')
    cols = {}
    for i, x in enumerate(line):
        cols[x] = i
    return cols


mztab_header = [
    "PSH", "PSM_ID", "accession", "unique", "database", "database_version",
    "search_engine", "search_engine_score[1]", "sequence", "modifications",
    "spectra_ref", "retention_time", "charge", "exp_mass_to_charge",
    "calc_mass_to_charge", "pre", "post", "start", "end"
]


def main():
    psm_file = sys.argv[1]
    output = sys.argv[2]
    fh = open(psm_file)
    out = open(output, 'w')
    out.write('\t'.join(mztab_header) + '\n')
    cols = get_cols(fh.readline())
    for line in fh:
        line = line.strip().split('\t')
        proteins = line[cols["protein"]].split(";")
        starts = line[cols["start"]].split(";")
        ends = line[cols["end"]].split(";")
        for accession, start, end in zip(proteins, starts, ends):
            if accession.startswith('sp') or accession.startswith('tr'):
                continue
            psm_id = line[cols['pride id']] + ':' + line[
                cols['sample']] + ':' + line[cols['feature_id']]

            sequence = line[cols['denovo_seq_nomod']]
            unique = '0'
            database = 'ENSEMBL'
            database_version = 'r94'
            search_engine = 'DeepNovoV2'
            search_engine_score = line[cols['predicted_score']]
            modifications = 'null'
            spectra_ref = line[cols['feature_id']]
            retention_time = 'null'
            charge = line[cols['precursor_charge']].split(".")[0]
            calc_mz = str(
                mass.calculate_mass(sequence=sequence, charge=int(charge)))
            charge += '+'
            exp_mass_to_charge = line[cols["precursor_mz"]]
            pre = 'null'
            post = 'null'
            out_list = [
                'PSM', psm_id, accession, unique, database, database_version,
                search_engine, search_engine_score, sequence, modifications,
                spectra_ref, retention_time, charge, exp_mass_to_charge,
                calc_mz, pre, post, start, end
            ]
            out.write('\t'.join(out_list) + '\n')
    out.close()
    fh.close()


if __name__ == "__main__":
    main()
