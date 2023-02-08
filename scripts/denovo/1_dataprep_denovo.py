import math
import re
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.style as style
import numpy as np
import pandas as pd

# in case you get a font warning please execute the following lines
# sudo apt-get install msttcorefonts -qq
# rm -rf ~/.cache/matplotlib/*

# figure parameters
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['axes.titlesize'] = 6
mpl.rcParams['axes.labelsize'] = 6
mpl.rcParams['xtick.labelsize'] = 6
mpl.rcParams['ytick.labelsize'] = 6
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['figure.autolayout'] = True  # tight_layout
mpl.rcParams['axes.spines.top'] = False
mpl.rcParams['axes.spines.right'] = False
mpl.rcParams['legend.frameon'] = False
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.sans-serif'] = "Arial"

denovo_file = sys.argv[1]
accuracy_file = sys.argv[2]

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
    'C': 103.00919,  # 4
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

COUNT_AA = {
    'A': [],
    'R': [],
    'N': [],
    'D': [],
    'C': [],
    'E': [],
    'Q': [],
    'G': [],
    'H': [],
    'I': [],
    'K': [],
    'M': [],
    'F': [],
    'P': [],
    'S': [],
    'T': [],
    'W': [],
    'Y': [],
    'V': []
}

LC_AA = {}


def splitAtUpperCase(text):
    result = ""
    for char in text:
        if char.isupper():
            result += " " + char
        else:
            result += char
    return result.split()


def AA_LC_CALCULATOR(decoder_input, output):
    """TODO(nh2tran): docstring."""
    if decoder_input == 'nan' or output == 'nan':
        return None
    decoder_input_len = len(splitAtUpperCase(decoder_input))
    output_len = len(splitAtUpperCase(output))
    decoder_input_mass = [mass_AA[x] for x in splitAtUpperCase(decoder_input)]
    decoder_input_mass_cum = np.cumsum(decoder_input_mass)
    output_mass = [mass_AA[x] for x in splitAtUpperCase(output)]
    output_mass_cum = np.cumsum(output_mass)
    i = 0
    j = 0
    while i < decoder_input_len and j < output_len:
        if abs(decoder_input_mass_cum[i] - output_mass_cum[j]) < 0.5:
            if decoder_input[i] == output[j]:
                COUNT_AA[decoder_input[i]].append(1)
            else:
                COUNT_AA[decoder_input[i]].append(0)
            i += 1
            j += 1
        elif decoder_input_mass_cum[i] < output_mass_cum[j]:
            COUNT_AA[decoder_input[i]].append(0)
            i += 1
        else:
            j += 1
    return None


def match_apply(row):
    if row.denovo_seq_nomod is not np.nan and row.db_search_seq is not np.nan:
        return test_AA_match_novor(row.denovo_seq_nomod, row.db_search_seq)
    else:
        return np.nan


def calc_lc(peptide):
    if isinstance(peptide, str):
        lc_list = []
        for aa in peptide:
            lc_list.append(str(round(LC_AA.get(aa, 0), 2) * 100))
        return ' '.join(lc_list)
    else:
        return 0


def mean(x):
    if len(x) == 0:
        return 0
    return sum(x) / len(x)


def calc_alc(peptide):
    if isinstance(peptide, str):
        lc_list = []
        for aa in peptide:
            lc_list.append(LC_AA.get(aa, 0) * 100)
        return mean(lc_list)
    else:
        return 0


def read_feature_accuracy(input_file, split_char):
    feature_list = []
    with open(input_file, 'r') as handle:
        header_line = handle.readline()
        for line in handle:
            line = re.split(split_char, line)
            feature = {}
            feature["feature_id"] = line[0]
            feature["feature_area"] = math.log10(float(line[1]) + 1e-5)
            feature["predicted_score"] = float(line[4])
            feature["recall_AA"] = float(line[5])
            feature["predicted_len"] = float(line[6])
            feature_list.append(feature)
    return feature_list


def find_score_cutoff(accuracy_file, accuracy_cutoff):
    print("".join(["="] * 80))  # section-separating line
    print("find_score_cutoff()")

    feature_list = read_feature_accuracy(accuracy_file, '\t|\r|\n')
    feature_list_sorted = sorted(
        feature_list, key=lambda k: k['predicted_score'], reverse=True)
    recall_cumsum = np.cumsum([f['recall_AA'] for f in feature_list_sorted])
    predicted_len_cumsum = np.cumsum(
        [f['predicted_len'] for f in feature_list_sorted])
    accuracy_cumsum = recall_cumsum / predicted_len_cumsum
    try:
        cutoff_index = np.flatnonzero(accuracy_cumsum < accuracy_cutoff)[0]
    except IndexError:
        cutoff_index = len(accuracy_cumsum) - 1
    cutoff_score = feature_list_sorted[cutoff_index]['predicted_score']
    while cutoff_score == -999:
        cutoff_index = cutoff_index - 1
        cutoff_score = feature_list_sorted[cutoff_index]['predicted_score']
    print('cutoff_score = ', cutoff_score)
    return cutoff_score


def convert_seq(seq):
    if isinstance(seq, str):
        if seq:
            return ''.join([char[0] for char in seq.split(',')])
    return seq


accuracy = pd.read_csv(accuracy_file, sep='\t', header=0, index_col=False)
denovo = pd.read_csv(denovo_file, sep='\t', header=0, index_col=False)

accuracy.predicted_sequence.fillna('', inplace=True)
denovo.predicted_sequence.fillna('', inplace=True)
denovo.predicted_score.fillna(-999, inplace=True)

score_cutoff = find_score_cutoff(accuracy_file, 0.90)

denovo['denovo_seq'] = denovo.predicted_sequence.str.replace(
    ',', '', regex=False)
denovo['denovo_seq_nomod'] = denovo.predicted_sequence.map(convert_seq)

cond1 = (pd.isnull(denovo.denovo_seq) == False)
cond2 = (denovo.denovo_seq_nomod.str.len() > 7)
cond3 = (denovo.denovo_seq_nomod.str.len() < 13)
cond4 = (denovo.predicted_score > score_cutoff)
denovo = denovo[cond1 & cond2 & cond3 & cond4]

denovo.to_csv('denovo_data_prep.tsv', sep='\t', header=True, index=False)

accuracy['denovo_seq'] = accuracy.predicted_sequence.str.replace(
    ',', '', regex=False)

accuracy['denovo_seq_nomod'] = accuracy.predicted_sequence.map(convert_seq)
accuracy['db_search_seq'] = accuracy.target_sequence.map(convert_seq)

accuracy['denovo_seq_nomod'] = accuracy['denovo_seq_nomod'].str.replace(
    'L', 'I')
accuracy['db_search_seq'] = accuracy['db_search_seq'].str.replace('L', 'I')
cond1 = pd.isnull(accuracy['target_sequence']) == False
cond2 = pd.isnull(accuracy['denovo_seq_nomod']) == False
srch_df = accuracy[(cond1) & (cond2)]

accuracy_list = srch_df['denovo_seq_nomod'].tolist()
srch_list = srch_df['db_search_seq'].tolist()

for denovo_pep, srch_pep in zip(accuracy_list, srch_list):
    AA_LC_CALCULATOR(denovo_pep, srch_pep)

for key, value in COUNT_AA.items():
    numerator = sum(value)
    denominator = len(value)
    if denominator == 0:
        LC_AA[key] = 0
    else:
        LC_AA[key] = sum(value) / len(value)

accuracy['LC'] = srch_df['denovo_seq_nomod'].map(calc_lc)
accuracy['ALC'] = accuracy['denovo_seq_nomod'].map(calc_alc)

style.use('seaborn-deep')
fig, ax = plt.subplots(1, 1, figsize=(5, 4))
pd.Series(LC_AA).plot(kind='bar', ax=ax)
plt.savefig('aa_LC.pdf', dpi=300, transparent=True)

# after applying the score cutoff
cond3 = accuracy.predicted_score >= score_cutoff
srch_df = accuracy[(cond1) & (cond2) & (cond3)]

accuracy_list = srch_df['denovo_seq_nomod'].tolist()
srch_list = srch_df['db_search_seq'].tolist()

for denovo_pep, srch_pep in zip(accuracy_list, srch_list):
    AA_LC_CALCULATOR(denovo_pep, srch_pep)

for key, value in COUNT_AA.items():
    numerator = sum(value)
    denominator = len(value)
    if denominator == 0:
        LC_AA[key] = 0
    else:
        LC_AA[key] = sum(value) / len(value)

style.use('seaborn-deep')
fig, ax = plt.subplots(1, 1, figsize=(5, 4))
pd.Series(LC_AA).plot(kind='bar', ax=ax)
plt.savefig('aa_LC_95pTP.pdf', dpi=300, transparent=True)
