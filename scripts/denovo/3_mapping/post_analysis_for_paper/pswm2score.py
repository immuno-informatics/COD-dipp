"""
    # Use example
    test = pswm["key"]
    test.max_score_()
    test.min_score_()
    test.calculate_score_("YANRNRFLY")
    test.calculate_scores_("YANRNRFLY")
    """
import argparse
import pickle as pkl
import re

import numpy as np
from Bio import SeqIO

letters = 'ACDEFGHIKLMNPQRSTVWYBZJO'


def max_score(pswm):
    m = pswm.length
    logodds = np.array([[pswm[letter][i] for letter in letters]
                        for i in range(m)], float)
    max_score = [logodds[i, :].max() for i in range(m)]
    return sum(max_score)


def min_score(pswm):
    m = pswm.length
    logodds = np.array([[pswm[letter][i] for letter in letters]
                        for i in range(m)], float)
    max_score = [logodds[i, :].min() for i in range(m)]
    return sum(max_score)


def calculate_scores(pswm, sequence):
    scores = [pswm[aa][i] for i, aa in enumerate(sequence)]
    return scores


def args():
    parser = argparse.ArgumentParser(
        description='calculates scores through pswm given a fasta input')
    requiredNamed = parser.add_argument_group('Required named arguments')
    requiredNamed.add_argument('-input_fasta',
                               nargs=1,
                               type=str,
                               help='Path to input peptides fasta',
                               required=True)
    requiredNamed.add_argument('-input_pswm',
                               nargs=1,
                               type=str,
                               help='Path to input pswm pickle object',
                               required=True)
    requiredNamed.add_argument('-output',
                               nargs=1,
                               type=str,
                               help='Path to output tsv',
                               required=True)
    optionalNamed = parser.add_argument_group('Optional named arguments')
    optionalNamed.add_argument('-sample_name',
                               nargs=1,
                               type=str,
                               default="all",
                               help='Use PSWMs from specific sample',
                               required=False)

    return parser


if __name__ == "__main__":
    arguments = args().parse_args()
    input_file = arguments.input_fasta[0]
    input_pswm = arguments.input_pswm[0]
    output_file = arguments.output[0]
    sample_name = arguments.sample_name[0]

    tmp_pswm = pkl.load(open(input_pswm, 'rb'))
    if sample_name != "all":
        pswm = {}
        for k, v in tmp_pswm.items():
            if sample_name in k:
                pswm[k] = v
    else:
        pswm = tmp_pswm

    if len(pswm) == 0:
        print("no PSMWs found in the pickle file")
        exit(1)

    peptides = SeqIO.index(input_file, format='fasta')
    output = open(output_file, 'w')

    header = 'peptide\tsample\tcluster\tlength\tpswm\tscores\t'
    header += 'score\tmax score\tmin score\tnormalized score\n'
    output.write(header)

    max_min_scores = {}
    keys = list(pswm.keys())
    lens = []

    for key, matrix in pswm.items():
        # key format = "{sample} cluster: {cluster} length: {length}"
        lens.append(int(key.split(': ')[-1]))
        max_min_scores[key] = [max_score(matrix), min_score(matrix)]
    lens = np.array(lens)
    max_len = max(lens)
    lin_len = min(lens)
    aas = '^[A|R|N|D|C|Q|E|G|H|I|L|K|M|F|P|S|T|W|Y|V]+$'
    for header, seq in peptides.items():
        seq = str(seq.seq)
        if not re.match(aas, seq):
            continue
        lenpep = len(seq)
        idx = np.where(lens == lenpep)[0]
        for ind in idx:
            try:
                scores = calculate_scores(pswm[keys[ind]], seq)
                score = sum(scores)
            except:
                print(f"{seq} ignored")
                continue
            scores = ' '.join([str(x) for x in scores])
            max_score_ = max_min_scores[keys[ind]][0]
            min_score_ = max_min_scores[keys[ind]][1]
            nscore = score / max_score_
            splitted = keys[ind].split(' ')
            sample = splitted[0]
            cluster = splitted[2]
            length = splitted[4]
            line = [
                seq, sample, cluster, length, keys[ind], scores, score,
                max_score_, min_score_, nscore
            ]
            line = [str(x) for x in line]
            output.write('\t'.join(line) + '\n')
    output.close()
