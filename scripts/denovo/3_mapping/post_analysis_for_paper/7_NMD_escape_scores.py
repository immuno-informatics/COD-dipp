import pandas as pd
import pybedtools
from numpy import dtype
import numpy as np
import matplotlib.pyplot as plt

file1 = '/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new/NMDescape/hg38_NMDetectiveA_Lindeboom_et_al.v2.gtf'
file2 = '/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new/NMDescape/hg38_NMDetectiveB_Lindeboom_et_al.v2.gtf'

mutation_list = [
    "COSV104417156", "COSV105148774", "COSV105262184", "COSV105262186",
    "COSV105274039", "COSV50463038", "COSV50466824", "COSV50468664",
    "COSV50477159", "COSV50510065", "COSV50510843", "COSV52591780",
    "COSV53683051", "COSV53915837", "COSV54156331", "COSV55237979",
    "COSV55266227", "COSV55825504", "COSV55826256", "COSV56467034",
    "COSV56468326", "COSV56470785", "COSV56472213", "COSV56473538",
    "COSV58346425", "COSV58388235", "COSV58584450", "COSV58727291",
    "COSV59364291", "COSV61086193", "COSV62425539", "COSV62426663",
    "COSV63114245", "COSV63114431", "COSV63479979", "COSV63631015",
    "COSV99064253", "COSV99064255", "COSV99064256", "COSV99180367",
    "COSV99802898", "COSV99875103"
]

cosmic_muts_file = "/net/archive/groups/plgg_iccvs/Resources/COSMIC/grch38_v92/CosmicMutantExport.tsv"

cols = ["GENOMIC_MUTATION_ID", "Mutation genome position"]

dtypes = {
    'Mutation genome position': dtype('O'),
}
print('reading redundant cosmic table entries')
df_cosmic = pd.read_csv(cosmic_muts_file,
                        sep="\t",
                        header=0,
                        usecols=cols,
                        dtype=dtypes)
df_cosmic = df_cosmic.drop_duplicates()
df_cosmic[df_cosmic['GENOMIC_MUTATION_ID'].isin(mutation_list)]
df_cosmic = df_cosmic[df_cosmic['GENOMIC_MUTATION_ID'].isin(mutation_list)]
df_cosmic[['chr', 'start', 'end'
           ]] = df_cosmic['Mutation genome position'].str.split(':|-',
                                                                expand=True)
df_cosmic['score'] = '.'
df_cosmic['strand'] = '.'

df1 = pd.read_csv(file1, sep='\t', header=None)
df2 = pd.read_csv(file2, sep='\t', header=None)
df2['score'] = '.'
df2['name'] = range(df2.shape[0])
df2[0] = df2[0].str.replace('^chr', '', regex=True)
with_scores = df2[8].str.split('; ', expand=True)
with_scores = with_scores[~with_scores.isnull().apply(any, axis=1)]

cols = [0, 3, 4, "name", "score", 6]
nmd_bed = pybedtools.BedTool.from_dataframe(df2.iloc[with_scores.index][cols])

cols = ['chr', 'start', 'end', 'GENOMIC_MUTATION_ID', 'score', 'strand']
fs_bed = pybedtools.BedTool.from_dataframe(df_cosmic[cols])

print("merging cosmic mutations with nmd scores")
inter = nmd_bed.intersect(fs_bed, loj=True).to_dataframe()
inter = inter[inter["thickStart"] != "."]
cols = ['chr', 'start', 'end', 'id', 'score', "strand"]
columns = [x + "_nmd" for x in cols]
columns.extend([x + "_fs" for x in cols])
inter.columns = columns


def get_scores(sdf, index, start, end):
    return [float(x) for x in sdf.loc[index, 4][10:-1].split(',')][start:end]


inter['range_start'] = inter.start_fs - inter.start_nmd
inter['range_end'] = inter.end_fs - inter.start_nmd
inter['nmd_score'] = inter.apply(lambda x: get_scores(
    with_scores, x.id_nmd, x.range_start, x.range_end + 1),
                                 axis=1)

fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
plot_series = inter.groupby('id_fs').apply(
    lambda x: x.nmd_score.map(np.mean).mean())
plot_series.plot(kind='hist', bins=10, edgecolor='white', ax=ax)
plt.savefig('nmd_efficacy_mut.pdf', transparent=True)

converter = {
    "COSV56467034": "ALFSRPTQPK",
    "COSV56468326": "ALFSRPTQPK",
    "COSV56470785": "ALFSRPTQPK",
    "COSV56472213": "ALFSRPTQPK",
    "COSV56473538": "ALFSRPTQPK",
    "COSV99802898": "ALFSRPTQPK",
    "COSV59364291": "ATHHQPWAQK",
    "COSV54156331": "AVMMNVRVAR",
    "COSV58346425": "DEFQLVSRY",
    "COSV105262184": "ESDLHPQKY",
    "COSV105262186": "ESDLHPQKY",
    "COSV61086193": "ESDLHPQKY",
    "COSV104417156": "GMPWHFMLK",
    "COSV105274039": "GMPWHFMLK",
    "COSV63479979": "GMPWHFMLK",
    "COSV105148774": "IVHTNFVEK",
    "COSV58727291": "KTWMMRKTW",
    "COSV63631015": "KVLDTPHHSK",
    "COSV55825504": "KVYCSFTRK",
    "COSV55826256": "KVYCSFTRK",
    "COSV99875103": "KVYCSFTRK",
    "COSV63114245": "MPKKSRISF",
    "COSV63114431": "MPKKSRISF",
    "COSV55266227": "QLWTVKLLK",
    "COSV58388235": "RLFNSVVPAYR",
    "COSV62425539": "RLWKHTLKY",
    "COSV62426663": "RLWKHTLKY",
    "COSV50463038": "RMSTAVYRW",
    "COSV50466824": "RMSTAVYRW",
    "COSV50468664": "RMSTAVYRW",
    "COSV50477159": "RMSTAVYRW",
    "COSV50510065": "RMSTAVYRW",
    "COSV50510843": "RMSTAVYRW",
    "COSV99064253": "RMSTAVYRW",
    "COSV99064255": "RMSTAVYRW",
    "COSV99064256": "RMSTAVYRW",
    "COSV52591780": "RPSPVRVAAL",
    "COSV99180367": "RPSPVRVAAL",
    "COSV58584450": "RSASYWPTK",
    "COSV53683051": "RTHQQAPLK",
    "COSV55237979": "RVQTSHTSWK",
    "COSV53915837": "SLSSWSLMK"
}

df_scores = plot_series.reset_index()
df_scores.columns = ['id_fs', 'nmd_score']
df_scores['peptide'] = df_scores.id_fs.map(converter)
df_scores.to_csv(
    '/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new/NMDescape/nmd_efficacy_scores.tsv',
    sep='\t',
    header=True,
    index=False)
df_scores_min = df_scores.groupby('peptide')['nmd_score'].min()
df_scores_min.to_csv(
    '/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new/NMDescape/nmd_efficacy_min_scores.tsv',
    sep='\t',
    header=True,
    index=False)
plot_series_pep = pd.cut(df_scores, bins=np.arange(0, 1.1, 0.1)).value_counts()

fig, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=True)
df_scores_min.plot(kind='hist',
                   bins=np.arange(0, 1.1, 0.1),
                   edgecolor='white',
                   ax=ax)
plt.savefig(
    '/net/archive/groups/plgg_iccvs/Datasets/proteomics_datasets/immunopeptidomics/PRIDE/meta_analysis/data_all/denovo_annotation/overlap_new/NMDescape/nmd_efficacy_pep.pdf',
    transparent=True)

df_
