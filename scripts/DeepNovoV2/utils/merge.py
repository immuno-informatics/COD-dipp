import csv
import os
import re
import sys

from pyteomics import mgf

import numpy as np

col_feature_id = 0
col_scan_list = 5
col_raw_sequence = 4

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


def cat_file_mgf(input_file_list, output_file):
    print("cat_file_mgf()")

    counter = 0
    with open(output_file, mode="w") as output:
        for index, input_file in enumerate(input_file_list):
            if 'spectrums.mgf' in input_file:
                continue
            print("input_file = ", os.path.join(input_file))
            with open(input_file) as input_handle:
                for line in input_handle:
                    if line.startswith('BEGIN IONS'):
                        line = input_handle.readline()
                        counter += 1
                        title = line.split(' ')[0] + '\n'
                        # SCANS=file:scan
                        scans = "SCANS=" + re.split('=|\n|"| ', line)[1] + ':'
                        scans = scans + title.rsplit('.')[-2]
                        rtinseconds = input_handle.readline()
                        pepmass = input_handle.readline().strip("\n").split(
                            ' ')[0] + "\n"
                        charge = input_handle.readline()
                        output.write("BEGIN IONS\n")
                        output.write(title)
                        output.write(pepmass)
                        output.write(charge)
                        output.write(scans + '\n')
                        output.write(rtinseconds)
                    else:
                        output.write(line)
    print("output_file = {0:s}".format(output_file))
    print("counter = {0:d}".format(counter))


def split_identified_and_unidentified_features(
        feature_file_name: str,
        output_identified_file_name=None,
        output_unidentified_file_name=None):
    if output_identified_file_name is None:
        output_identified_file_name = feature_file_name + '.identified'
    if output_unidentified_file_name is None:
        output_unidentified_file_name = feature_file_name + '.unidentified'
    id_handle = open(output_identified_file_name, 'w')
    unid_handle = open(output_unidentified_file_name, 'w')

    with open(feature_file_name, 'r') as f:
        line = f.readline()
        # write header
        id_handle.write(line)
        unid_handle.write(line)
        for line in f:
            seq = line.split(',')[col_raw_sequence]
            if seq:
                id_handle.write(line)
            else:
                unid_handle.write(line)


def partition_feature_file_nodup(input_feature_file, prob):
    print("partition_feature_file_nodup()")

    print("input_feature_file = ", os.path.join(input_feature_file))
    print("prob = ", prob)

    output_file_train = input_feature_file + ".train" + ".nodup"
    output_file_valid = input_feature_file + ".valid" + ".nodup"
    output_file_test = input_feature_file + ".test" + ".nodup"

    peptide_train_list = []
    peptide_valid_list = []
    peptide_test_list = []

    with open(input_feature_file, mode="r") as input_handle:
        with open(output_file_train, mode="w") as output_handle_train:
            with open(output_file_valid, mode="w") as output_handle_valid:
                with open(output_file_test, mode="w") as output_handle_test:
                    counter = 0
                    counter_train = 0
                    counter_valid = 0
                    counter_test = 0
                    counter_unique = 0
                    # header line
                    line = input_handle.readline()
                    output_handle_train.write(line)
                    output_handle_valid.write(line)
                    output_handle_test.write(line)
                    # first feature
                    line = input_handle.readline()
                    while line and counter_unique < 4000:
                        counter += 1
                        # check if the peptide already exists in any of the three lists
                        # if yes, this new feature will be assigned to that list
                        peptide = re.split(',|\r|\n', line)[col_raw_sequence]
                        if (peptide in peptide_train_list):
                            output_handle = output_handle_train
                            counter_train += 1
                        elif (peptide in peptide_valid_list):
                            output_handle = output_handle_valid
                            counter_valid += 1
                        elif (peptide in peptide_test_list):
                            output_handle = output_handle_test
                            counter_test += 1
                        # if not, this new peptide and its spectrum will be randomly assigned
                        else:
                            counter_unique += 1
                            set_num = np.random.choice(a=3, size=1, p=prob)
                            if set_num == 0:
                                peptide_train_list.append(peptide)
                                output_handle = output_handle_train
                                counter_train += 1
                            elif set_num == 1:
                                peptide_valid_list.append(peptide)
                                output_handle = output_handle_valid
                                counter_valid += 1
                            else:
                                peptide_test_list.append(peptide)
                                output_handle = output_handle_test
                                counter_test += 1
                        output_handle.write(line)
                        line = input_handle.readline()

    input_handle.close()
    output_handle_train.close()
    output_handle_valid.close()
    output_handle_test.close()

    print("counter = {0:d}".format(counter))
    print("counter_train = {0:d}".format(counter_train))
    print("counter_valid = {0:d}".format(counter_valid))
    print("counter_test = {0:d}".format(counter_test))
    print("counter_unique = {0:d}".format(counter_unique))


def partition_mgf_file_nodup(input_mgf_file, input_feature_file_train,
                             input_feature_file_valid,
                             input_feature_file_test):
    feature_list = [
        input_feature_file_train, input_feature_file_valid,
        input_feature_file_test
    ]
    dir_name = os.path.dirname(input_feature_file_train)

    mgf_list = [
        os.path.join(dir_name, x)
        for x in ['training.mgf', 'validation.mgf', 'testing.mgf']
    ]

    mf = mgf.read(input_mgf_file, use_index=True)

    for feature_file, mgf_file in zip(feature_list, mgf_list):
        spectra_dict = {}
        with open(feature_file) as ff:
            header = ff.readline()
            for line in ff:
                line = line.strip('\n').split(',')
                spectrum_id = line[-1]
                spectra_dict[spectrum_id] = True
        output = open(mgf_file, 'w')
        with open(input_mgf_file) as input_handle:
            line = input_handle.readline()
            while line:
                if line.startswith('BEGIN IONS'):
                    title = re.split('=|\n', input_handle.readline())[1]
                    if spectra_dict.get(title, False):
                        title = 'TITLE=' + title + '\n'
                        # SCANS=file:scan
                        scans = input_handle.readline()
                        rtinseconds = input_handle.readline()
                        pepmass = input_handle.readline()
                        charge = input_handle.readline()
                        output.write("BEGIN IONS\n")
                        output.write(title)
                        output.write(pepmass)
                        output.write(charge)
                        output.write(scans)
                        output.write(rtinseconds)
                        line = charge
                    else:
                        while line:
                            line = input_handle.readline()
                            if line.startswith('BEGIN ION'):
                                break
                else:
                    line = input_handle.readline()
                    output.write(line)

        output.close()


def feature_extract(search_engine, post_processing, input_psm_file, mgf_file,
                    output_file):
    col_spectrum_title = 'Spectrum'
    col_peptide = 'Peptide'
    col_obs_mod = 'Observed Modifications'
    col_ass_mod_tpp = 'Assigned Modifications'
    col_ass_mod_scavager = 'modifications'
    col_exp_mass = 'Experimental Mass'
    col_theo_mass = 'Peptide Mass'
    col_ppm_diff = 'massdiff_ppm'

    mf = mgf.read(mgf_file)
    header = [
        "spec_group_id", "m/z", "z", "rt_mean", "seq", "scans", "profile",
        "feature area", 'title'
    ]
    if post_processing == 'scavager':
        col_spectrum_title = col_spectrum_title.lower()
        col_peptide = col_peptide.lower()
    with open(output_file, 'w') as output_handle:
        output_handle.write(','.join(header))
        output_handle.write('\n')
        with open(input_psm_file) as input_handle:
            seq_dict = {}
            header = input_handle.readline().strip('\n').split('\t')
            cols = {}
            for index, col in enumerate(header):
                cols[col] = index
            for line in input_handle:
                line = line.strip('\n').split('\t')
                try:
                    [mass_AA[x] for x in line[cols[col_peptide]]]
                except KeyError:
                    continue
                if post_processing == 'tpp':
                    observed_mass = float(line[cols[col_exp_mass]])
                    theoretical_mass = float(line[cols[col_theo_mass]])
                    ppm = (observed_mass -
                           theoretical_mass) / theoretical_mass * 1e6
                    if line[cols[col_ass_mod_tpp]] != '':
                        continue
                elif post_processing == 'scavager':
                    ppm = float(line[cols[col_ppm_diff]])
                    try:
                        if line[cols[
                                col_ass_mod_scavager]] != '[]' and search_engine == 'msfragger':
                            continue
                        if line[cols[
                                col_ass_mod_scavager]] != '' and search_engine == 'msgfplus':
                            continue
                    except KeyError:
                        pass
                else:
                    raise ValueError(
                        f"{post_processing} post processing unexpected")
                if abs(ppm) <= 10:
                    spec_title = line[cols[col_spectrum_title]]
                    spec_title = re.sub('\\.0+', '\\.', spec_title)
                    seq = line[cols[col_peptide]]
                    seq_dict[spec_title] = seq
        for spectrum in mf:
            spec_title = spectrum['params']['title']
            seq = seq_dict.get(spec_title, '')
            try:
                row = [
                    spectrum['params']['scans'],
                    spectrum['params']['pepmass'][0],
                    spectrum['params']['charge'][0],
                    spectrum['params']['rtinseconds'], seq,
                    spectrum['params']['scans'], "0.0:1.0", "1.0", spec_title
                ]
            except:
                print("{} ignored".format(spec_title))
                continue
            output_handle.write(','.join(
                [str(item).strip('+') for item in row]))
            output_handle.write('\n')
    mf.close()
