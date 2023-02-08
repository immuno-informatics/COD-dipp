import csv
import glob
import os
import sys

import numpy as np
from utils.merge import (cat_file_mgf, feature_extract,
                         partition_feature_file_nodup,
                         partition_mgf_file_nodup,
                         split_identified_and_unidentified_features)

mass_H = 1.0078
mass_H2O = 18.0106
mass_NH3 = 17.0265
mass_N_terminus = 1.0078
mass_C_terminus = 17.0027
mass_CO = 27.9949

mass_AA = {
    '_PAD': 0.0,
    '_GO': mass_N_terminus - mass_H,
    '_EOS': mass_C_terminus + mass_H,
    'A': 71.03711,  # 0
    'R': 156.10111,  # 1
    'N': 114.04293,  # 2
    'D': 115.02694,  # 3
    'C': 103.00918,
    'E': 129.04259,  # 5
    'Q': 128.05858,  # 6
    'G': 57.02146,  # 7
    'H': 137.05891,  # 8
    'I': 113.08406,  # 9
    'L': 113.08406,  # 10
    'K': 128.09496,  # 11
    'M': 131.04049,  # 12
    'F': 147.06841,  # 13
    'P': 97.05276,  # 14
    'S': 87.03203,  # 15
    'T': 101.04768,  # 16
    'W': 186.07931,  # 17
    'Y': 163.06333,  # 18
    'V': 99.06841,  # 19
}


def parse_raw_sequence(raw_sequence: str):
    raw_sequence_len = len(raw_sequence)
    peptide = []
    index = 0
    while index < raw_sequence_len:
        if raw_sequence[index] == "(":
            if peptide[-1] == "C" and raw_sequence[index:index +
                                                   8] == "(+57.02)":
                peptide[-1] = "C(Carbamidomethylation)"
                index += 8
            elif peptide[-1] == 'M' and raw_sequence[index:index +
                                                     8] == "(+15.99)":
                peptide[-1] = 'M(Oxidation)'
                index += 8
            elif peptide[-1] == 'N' and raw_sequence[index:index +
                                                     6] == "(+.98)":
                peptide[-1] = 'N(Deamidation)'
                index += 6
            elif peptide[-1] == 'Q' and raw_sequence[index:index +
                                                     6] == "(+.98)":
                peptide[-1] = 'Q(Deamidation)'
                index += 6
            else:  # unknown modification
                logger.warning(f"unknown modification in seq {raw_sequence}")
                return False, peptide
        else:
            peptide.append(raw_sequence[index])
            index += 1

    return True, peptide


def compute_neutral_peptide_mass(peptide: list):
    peptide_neutral_mass = mass_N_terminus + mass_C_terminus
    for aa in peptide:
        peptide_neutral_mass += mass_AA[aa]
    return peptide_neutral_mass


def feature_file_mass_correction(feature_filename: str):
    """
    read feature file, find out mass shift then correct
    :param feature_filename:
    :return:
    """
    output_feature_filename = feature_filename + '.mass_corrected'
    ppm_shift = []
    with open(feature_filename, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        header = next(reader)
        seq_index = header.index("seq")
        mz_index = header.index("m/z")
        z_index = header.index("z")
        for line in reader:
            mz = float(line[mz_index])
            z = float(line[z_index])
            observed_mass = mz * z - z * mass_H
            if not line[seq_index]:
                continue
            okay, peptide = parse_raw_sequence(line[seq_index])
            if not okay:
                # unknown mods
                continue
            theoretical_mass = compute_neutral_peptide_mass(peptide)
            ppm = (observed_mass - theoretical_mass) / theoretical_mass * 1e6
            ppm_shift.append(ppm)
    if len(ppm_shift) < 100:
        raise ValueError("too less identified feature for mass correction")
    ppm_shift = np.median(ppm_shift)
    print(f"ppm shift: {ppm_shift}")
    with open(feature_filename, 'r') as fr:
        with open(output_feature_filename, 'w') as fw:
            reader = csv.reader(fr, delimiter=',')
            writer = csv.writer(fw, delimiter=',')
            writer.writerow(next(reader))
            for line in reader:
                mz = float(line[mz_index])
                mz = mz * (1 - ppm_shift * 1e-6)
                line[mz_index] = "{}".format(mz)
                writer.writerow(line)


if __name__ == '__main__':
    prefix_folder_name = sys.argv[1]
    psm_file = sys.argv[2]
    post_processing = sys.argv[3]  # 'scavager' or 'tpp'
    search_engine = sys.argv[4]  # 'msgfplus' or 'msfragger'

    input_mgf_file_list = glob.glob(os.path.join(prefix_folder_name, '*.mgf'))
    feature_output_file = os.path.join(prefix_folder_name, 'export.csv')

    mgf_output_file = os.path.join(prefix_folder_name, 'spectrums.mgf')

    # concatenate different fractions
    cat_file_mgf(input_mgf_file_list, mgf_output_file)

    # extract features from the db search results
    feature_extract(search_engine, post_processing, psm_file, mgf_output_file,
                    feature_output_file)

    # mass correction
    feature_file_mass_correction(feature_output_file)

    split_identified_and_unidentified_features(feature_output_file +
                                               '.mass_corrected')

    # split identified features for training, validation and testing
    partition_feature_file_nodup(
        feature_output_file + '.mass_corrected.identified', [0.8, 0.1, 0.1])

    # split concatenated mgf file for training, validation and testing
    feature_file_train = feature_output_file + '.mass_corrected.identified.train.nodup'
    feature_file_valid = feature_output_file + '.mass_corrected.identified.valid.nodup'
    feature_file_test = feature_output_file + '.mass_corrected.identified.test.nodup'

    # partition_mgf_file_nodup(mgf_output_file, feature_file_train,
    #                          feature_file_valid, feature_file_test)
