from __future__ import division
import argparse
# from . import main, utils
import pkg_resources
import pandas as pd
import ast
import subprocess
import numpy as np

def run():
    parser = argparse.ArgumentParser(
        description='run diffacto for scavager results',
        epilog='''

    Example usage
    -------------
    $ scav2diffacto -S1 sample1_1_proteins.tsv sample1_n_proteins.tsv -S2 sample2_1_proteins.tsv sample2_n_proteins.tsv
    -------------
    ''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-dif', help='path to Diffacto', required=True)
    parser.add_argument('-S1', nargs='+', help='input files for S1 sample', required=True)
    parser.add_argument('-S2', nargs='+', help='input files for S2 sample', required=True)
    parser.add_argument('-S3', nargs='+', help='input files for S3 sample')
    parser.add_argument('-S4', nargs='+', help='input files for S4 sample')
    parser.add_argument('-peptides', help='name of output peptides file', default='peptides.txt')
    parser.add_argument('-samples', help='name of output samples file', default='sample.txt')
    parser.add_argument('-out', help='name of diffacto output file', default='diffacto_out.txt')
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


    for sample_num in ['S1', 'S2', 'S3', 'S4']:
        if args[sample_num]:
            for z in args[sample_num]:
                df0 = pd.read_table(z.replace('_proteins.tsv', '_peptides.tsv'))
                allowed_peptides.update(df0['peptide'])


    for sample_num in ['S1', 'S2', 'S3', 'S4']:
        if args[sample_num]:
            for z in args[sample_num]:
                label = z.replace('_proteins.tsv', '')
                df1 = pd.read_table(z.replace('_proteins.tsv', '_PSMs_full.tsv'))
                # df1 = pd.read_table(z.replace('_proteins.tsv', '_PSMs.tsv'))
                df1 = df1[df1['peptide'].apply(lambda z: z in allowed_peptides)]
                # print(df1.shape)
                # print(z.replace('_proteins.tsv', '_PSMs_full.tsv'))
                # print(df1.columns)
                df1['peptide'] = df1.apply(lambda z: z['peptide'] + str(z['assumed_charge']), axis=1)
                df1 = df1.sort_values('MS1Intensity', ascending=False).drop_duplicates(['peptide'])
                df1['peptide'] = df1['peptide'].apply(lambda z: z[:-1])
                df1['MS1Intensity'] = df1.groupby('peptide')['MS1Intensity'].transform(sum)
                df1 = df1.sort_values('q', ascending=True).drop_duplicates(['peptide'])
                df1[label] = df1['MS1Intensity']
                df1[label] = df1[label].replace([0, 0.0], np.nan)
                df1['protein'] = df1['protein'].apply(lambda z: ';'.join([u for u in ast.literal_eval(z) if u in allowed_prots]))
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
    df_final.to_csv(args['peptides'], sep=',')

    out = open(args['samples'], 'w')
    for sample_num in ['S1', 'S2', 'S3', 'S4']:
        if args[sample_num]:
            for z in args[sample_num]:
                label = z.replace('_proteins.tsv', '')
                out.write(label + '\t' + sample_num + '\n')
    out.close()

    subprocess.call(['python3', args['dif'], '-i', args['peptides'], '-samples', args['samples'], '-out',\
     args['out'], '-normalize', args['norm'], '-impute_threshold', args['impute_threshold'], '-min_samples', args['min_samples']])



if __name__ == '__main__':
    run()
