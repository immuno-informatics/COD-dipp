from __future__ import division
import argparse
# from . import main, utils
import pkg_resources
import pandas as pd
import ast
import subprocess
import numpy as np
from itertools import chain
from pyteomics import fasta

def run():
    parser = argparse.ArgumentParser(
        description='run diffacto for scavager results',
        epilog='''

    Example usage
    -------------
    $ calcnsaf -S1 sample1_1_proteins.tsv sample1_n_proteins.tsv -S2 sample2_1_proteins.tsv sample2_n_proteins.tsv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-S1', nargs='+', help='input files for S1 sample', required=True)
    parser.add_argument('-S2', nargs='+', help='input files for S2 sample', required=True)
    parser.add_argument('-S3', nargs='+', help='input files for S3 sample')
    parser.add_argument('-S4', nargs='+', help='input files for S4 sample')
    parser.add_argument('-db', help='path to fasta file', required=True)
    parser.add_argument('-peptides', help='name of output peptides file', default='peptides.txt')
    parser.add_argument('-samples', help='name of output samples file', default='sample.txt')
    parser.add_argument('-out', help='name of nsaf output file', default='nsaf_out.txt')
    parser.add_argument('-norm', help='normalization method. Can be average, median, GMM or None', default='None')
    parser.add_argument('-impute_threshold', help='impute_threshold for missing values fraction', default='0.25')
    parser.add_argument('-min_samples', help='minimum number of samples for peptide usage', default='3')
    parser.add_argument('-version', action='version', version='%s' % (pkg_resources.require("scavager")[0], ))
    args = vars(parser.parse_args())


    df_final = False

    allowed_prots = set()
    allowed_peptides = set()

    for sample_num in ['S1', 'S2', 'S3', 'S4']:
        if args[sample_num]:
            for z in args[sample_num]:
                df0 = pd.read_table(z)
                allowed_prots.update(df0['dbname'])


    # for sample_num in ['S1', 'S2', 'S3', 'S4']:
    #     if args[sample_num]:
    #         for z in args[sample_num]:
    #             df0 = pd.read_table(z.replace('_proteins.tsv', '_peptides.tsv'))
    #             allowed_peptides.update(df0['peptide'])


    for sample_num in ['S1', 'S2', 'S3', 'S4']:
        if args[sample_num]:
            for z in args[sample_num]:
                label = z.replace('_proteins.csv', '')
                df1 = pd.read_table(z.replace('_proteins.csv', '_PFMs.csv'))
                # df1 = pd.read_table(z.replace('_proteins.tsv', '_PSMs.tsv'))
                df1['peptide'] = df1['sequence']
                # df1 = df1[df1['peptide'].apply(lambda z: z in allowed_peptides)]
                # print(df1.shape)
                # print(z.replace('_proteins.tsv', '_PSMs_full.tsv'))
                # print(df1.columns)
                # df1['peptide'] = df1.apply(lambda z: z['peptide'] + str(z['assumed_charge']), axis=1)
                # df1 = df1.sort_values('MS1Intensity', ascending=False).drop_duplicates(['peptide'])
                # df1['peptide'] = df1['peptide'].apply(lambda z: z[:-1])
                df1['count'] = df1.groupby('peptide')['peptide'].transform('count')
                # print(df1['count'])
                df1 = df1.drop_duplicates(['peptide'])
                df1[label] = df1['count']
                # df1[label] = df1[label].replace([0, 0.0], np.nan)
                df1['protein'] = df1['proteins'].apply(lambda z: ';'.join([u for u in z.split(';') if u in allowed_prots]))
                # df1['protein'] = df1['protein'].apply(lambda z: ';'.join([u for u in ast.literal_eval(z) if u in allowed_prots]))
                df1 = df1[df1['protein'].apply(lambda z: z != '')]
                df1 = df1[['peptide', 'protein', label]]
                if df_final is False:
                    df_final = df1
                else:
                    df_final = df_final.reset_index().merge(df1.reset_index(), on='peptide', how='outer')#.set_index('peptide')
                    # df_final = df_final.merge(df1, on='peptide', how='outer')
                    df_final.protein_x.fillna(value=df_final.protein_y, inplace=True)
                    df_final['protein'] = df_final['protein_x']
                    df_final = df_final.drop(columns=['protein_x', 'protein_y', 'index_x', 'index_y'])


    print(df_final.columns)
    df_final = df_final.set_index('peptide')
    df_final['proteins'] = df_final['protein']
    df_final = df_final.drop(columns=['protein'])
    cols = df_final.columns.tolist()
    cols.remove('proteins')
    cols.insert(0, 'proteins')
    df_final = df_final[cols]
    df_final.fillna(value='')

    cols = df_final.columns.difference(['proteins'])
    genres = df_final['proteins'].str.split(';')
    df_final =  (df_final.loc[df_final.index.repeat(genres.str.len()), cols]
         .assign(dbname=list(chain.from_iterable(genres.tolist()))))
    df_final = df_final.groupby('dbname').sum()
    df_final.reset_index(level=0, inplace=True)

    protsL = {}
    for p in fasta.read(args['db']):
        dbn = p[0].split()[0]
        protsL[dbn] = len(p[1])

    df_final['Length'] = df_final['dbname'].apply(lambda z: protsL[z])
    for cc in df_final.columns:
        if cc not in ['dbname', 'Length']:
            df_final[cc] = df_final[cc] / df_final['Length']
    for cc in df_final.columns:
        if cc not in ['dbname', 'Length']:
            df_final[cc] = df_final[cc] / df_final[cc].sum()
            df_final[cc] = df_final[cc].replace(0, np.nan)
            min_val = np.nanmin(df_final[cc].values)
            df_final[cc] = df_final[cc].replace(np.nan, min_val)
    df_final.drop(columns=['Length', ], inplace=True)
    df_final.to_csv(args['out'], sep='\t', index=False)

    print(df_final[:5])

    # out = open(args['samples'], 'w')
    # for sample_num in ['S1', 'S2', 'S3', 'S4']:
    #     if args[sample_num]:
    #         for z in args[sample_num]:
    #             label = z.replace('_proteins.tsv', '')
    #             out.write(label + '\t' + sample_num + '\n')
    # out.close()


if __name__ == '__main__':
    run()
