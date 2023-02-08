from Bio import SeqIO
import re
import pandas as pd
import os
root = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all"
file_ann = os.path.join(root, "tables/HLA-I_tables/samples_info.tsv")

db = '/net/archive/groups/plgg_iccvs/Resources/COSMIC/grch38_v92/protein_database/cosmic_proteins_v94.fasta'
file_input_cryp = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/denovo_nonexons_spectra_3m.tsv"
file_out = "/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new/stop_loss_peptides.tsv"

df_cryptic = pd.read_csv(file_input_cryp, sep='\t', header=0)
df_ann = pd.read_csv(file_ann, sep="\t")
df_ann = df_ann[df_ann["HLA Class"] == 1]

# search for stop loss mutations
pattern = re.compile(r".+p\.\*[0-9]+[^=]+$")
db_dict = {}

for seq in SeqIO.parse(db, format='fasta'):
    header = seq.id
    if pattern.match(header):
        seq = str(seq.seq)
        db_dict[header] = seq

peptides = list(df_cryptic.denovo_seq_nomod.unique())

found = {}

for k, v in db_dict.items():
    for pep in peptides:
        if pep.replace('I', 'L') in v.replace('I', 'L'):
            found[pep] = k

df = pd.Series(found).reset_index()
df.columns = ['Cryptic peptide', 'mutation']

df_count = df_cryptic[df_cryptic['denovo_seq_nomod'].isin(
    df['Cryptic peptide'])].groupby(
        ['denovo_seq_nomod',
         'sample']).apply(lambda x: x.shape[0]).reset_index().reset_index()
df_count.drop('index', axis=1, inplace=True)
df_count.columns = ['Cryptic peptide', 'Sample Name', '# PSM']
df_count2 = df.merge(df_count)

df_count3 = df_count2.merge(df_ann)

df_count3.to_csv(file_out, sep='\t', header=True, index=False)
