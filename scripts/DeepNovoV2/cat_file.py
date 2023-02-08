import csv
import glob
import os
import sys

import numpy as np

import deepnovo_config
from data_reader import parse_raw_sequence
from utils.merge import (cat_file_mgf, feature_extract,
                         partition_feature_file_nodup,
                         partition_mgf_file_nodup,
                         split_identified_and_unidentified_features)


def compute_neutral_peptide_mass(peptide: list):
    peptide_neutral_mass = deepnovo_config.mass_N_terminus + deepnovo_config.mass_C_terminus
    for aa in peptide:
        peptide_neutral_mass += deepnovo_config.mass_AA[aa]
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
            observed_mass = mz * z - z * deepnovo_config.mass_H
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

    input_mgf_file_list = glob.glob(os.path.join(prefix_folder_name, '*.mgf'))
    feature_output_file = os.path.join(prefix_folder_name, 'export.csv')

    mgf_output_file = os.path.join(prefix_folder_name, 'spectrum.mgf')

    # concatenate different fractions
    cat_file_mgf(input_mgf_file_list, mgf_output_file)

    # extract features from the db search results
    feature_extract(psm_file, mgf_output_file, feature_output_file)

    # mass correction
    feature_file_mass_correction(feature_output_file)

    # split_identified_and_unidentified_features(feature_output_file +
    #                                            '.mass_corrected')
    split_identified_and_unidentified_features(feature_output_file +
                                               '.mass_corrected')

    # split identified features for training, validation and testing
    partition_feature_file_nodup(
        feature_output_file + '.mass_corrected.identified', [0.8, 0.1, 0.1])

    # split concatenated mgf file for training, validation and testing
    suffix = '.mass_corrected.identified.train.nodup'
    feature_file_train = feature_output_file + suffix
    feature_file_valid = feature_output_file + suffix
    feature_file_test = feature_output_file + suffix

    partition_mgf_file_nodup(feature_file_train, feature_file_valid,
                             feature_file_test, mgf_output_file)
