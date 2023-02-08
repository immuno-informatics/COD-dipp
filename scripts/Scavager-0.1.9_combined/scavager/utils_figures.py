from __future__ import division
from scipy.stats import scoreatpercentile
from pyteomics import mass
from collections import Counter
import os.path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
try:
    import seaborn
    seaborn.set(rc={'axes.facecolor':'#ffffff'})
    seaborn.set_style('whitegrid')
except ImportError:
    pass
import re

redcolor = '#FC6264'
bluecolor = '#70aed1'
greencolor = '#8AA413'

def get_basic_distributions(df):
    mz_array = ((df['calc_neutral_pep_mass'] + df['assumed_charge'] * 1.007276 ) / df['assumed_charge']).values
    rt_exp_array = (df['RT exp']).values
    lengths_array = df['length'].values
    return mz_array, rt_exp_array, lengths_array

def get_descriptor_array(df, df_f, dname):
    array_t = df[~df['decoy']][dname].values
    array_d = df[df['decoy']][dname].values
    array_v = df_f[dname].values
    return array_t, array_d, array_v

def plot_hist_basic(array_all, array_valid, fig, subplot_max_x, subplot_i, xlabel, ylabel='# of identifications', bin_size_one=False):
    separate_figures = isinstance(fig, str)
    if separate_figures:
        plt.figure()
    else:
        fig.add_subplot(subplot_max_x, 3, subplot_i)
    cbins = get_bins((array_all, array_valid), bin_size_one)
    plt.hist(array_all, bins=cbins, color=redcolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.hist(array_valid, bins=cbins, color=greencolor, alpha=0.8, edgecolor='#EEEEEE')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if separate_figures:
        plt.savefig(outpath(fig, xlabel, '.png'))

def plot_basic_figures(df, df_f, fig, subplot_max_x, subplot_start, idtype):
    mz_array, rt_exp_array, lengths_array = get_basic_distributions(df)
    mz_array_valid, rt_exp_array_valid, lengths_array_valid = get_basic_distributions(df_f)

    plot_hist_basic(mz_array, mz_array_valid, fig, subplot_max_x, subplot_i=subplot_start, \
                     xlabel='%s, precursor m/z' % (idtype, ))
    subplot_start += 1
    plot_hist_basic(rt_exp_array, rt_exp_array_valid, fig, subplot_max_x, subplot_i=subplot_start, \
                     xlabel='%s, RT experimental' % (idtype, ))
    subplot_start += 1
    plot_hist_basic(lengths_array, lengths_array_valid, fig, subplot_max_x, subplot_i=subplot_start, \
                     xlabel='%s, peptide length' % (idtype, ), bin_size_one=True)

def plot_protein_figures(df, df_f, fig, subplot_max_x, subplot_start):
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='LOG10_NSAF'), fig, subplot_max_x, subplot_start, xlabel='LOG10(NSAF)')
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='sq'), fig, subplot_max_x, subplot_start+1, xlabel='sequence coverage')

def plot_hist_descriptor(inarrays, fig, subplot_max_x, subplot_i, xlabel, ylabel='# of identifications'):
    separate_figures = isinstance(fig, str)
    if separate_figures:
        plt.figure()
    else:
        fig.add_subplot(subplot_max_x, 3, subplot_i)
    array_t, array_d, array_v = inarrays
    cbins, width = get_bins_for_descriptors((array_t, array_d, array_v))
    H1, _ = np.histogram(array_d, bins=cbins)
    H2, _ = np.histogram(array_t, bins=cbins)
    H3, _ = np.histogram(array_v, bins=cbins)
    plt.bar(cbins[:-1], H1, width, align='center',color=redcolor, alpha=0.4, edgecolor='#EEEEEE')
    plt.bar(cbins[:-1], H2, width, align='center',color=bluecolor, alpha=0.4, edgecolor='#EEEEEE')
    plt.bar(cbins[:-1], H3, width, align='center',color=greencolor, alpha=1, edgecolor='#EEEEEE')
    cbins = np.append(cbins[0], cbins)
    H1 = np.append(np.append(0, H1), 0)
    H2 = np.append(np.append(0, H2), 0)
    cbins -= width / 2
    plt.step(cbins, H2, where='post', color=bluecolor, alpha=0.8)
    plt.step(cbins, H1, where='post', color=redcolor, alpha=0.8)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    if 'mass shift' in xlabel:
        plt.xlim(-1.5, 1.5)
        plt.xticks([-1, 0, 1], ['NA', 'unmodified', 'modified'])
    elif width == 1.0:
        plt.xticks(np.arange(int(cbins[0]), cbins[-1], 1))
        plt.gcf().canvas.draw()
    if separate_figures:
        plt.savefig(outpath(fig, xlabel, '.png'))

def plot_legend(fig, subplot_max_x, subplot_start):
    ax = fig.add_subplot(subplot_max_x, 3, subplot_start)
    legend_elements = [Patch(facecolor=greencolor, edgecolor='r',
                            label='Positive IDs'),
                    Patch(facecolor=bluecolor, edgecolor='r',
                            label='Targets'),
                    Patch(facecolor=redcolor, edgecolor='r',
                            label='Decoys')]
    ax.legend(handles=legend_elements, loc='center', prop={'size': 24})
    ax.set_axis_off()

def plot_aa_stats(df_f, df_proteins_f, fig, subplot_max_x, subplot_i):
    ax = fig.add_subplot(subplot_max_x, 3, subplot_i)

    # Generate list of 20 standart amino acids
    std_aa_list = list(mass.std_aa_mass.keys())
    std_aa_list.remove('O')
    std_aa_list.remove('U')

    # Count identified amino acids
    aa_exp = Counter()
    for pep in set(df_f['peptide']):
        for aa in pep:
            aa_exp[aa] += 1

    # Count theoretical amino acids
    aa_theor = Counter()
    for prot_seq in set(df_proteins_f['sequence'].values):
        for aa in prot_seq:
            aa_theor[aa] += 1

    aa_exp_sum = sum(aa_exp.values())
    aa_theor_sum = sum(aa_theor.values())
    lbls, vals = [], []
    for aa in sorted(std_aa_list):
        lbls.append(aa)
        vals.append((aa_exp.get(aa, 0)/aa_exp_sum)/(aa_theor.get(aa, 0)/aa_theor_sum))
    std_val = np.std(vals)
    clrs = [greencolor if abs(x-1)<=2*std_val else redcolor for x in vals]
    ax.bar(range(len(vals)), vals, color=clrs)
    ax.set_xticks(range(len(lbls)))
    ax.set_xticklabels(lbls)
    ax.hlines(1.0, range(len(vals))[0]-1, range(len(vals))[-1]+1)
    ax.set_ylabel('amino acid ID rate')


def calc_max_x_value(df, df_proteins):
    cnt = 7 # number of basic figures
    peptide_columns = set(df.columns)
    features_list = ['massdiff_ppm', 'RT diff', 'fragmentMT', 'num_missed_cleavages', 'assumed_charge', 'log_score', 'ISOWIDTHDIFF', 'MS1Intensity']
    for feature in features_list:
        if feature in peptide_columns:
            cnt += 1
    for feature in peptide_columns:
        if feature.startswith('mass shift'):
            cnt += 1
    if len(set(df['massdiff_int'])) > 1:
        cnt += 1
    if 'LOG10_NSAF' in df_proteins.columns:
        cnt += 3 # add for NSAF, sequence coverage and aa_stats
    return cnt // 3 + (1 if (cnt % 3) else 0)


def plot_descriptors_figures(df, df_f, fig, subplot_max_x, subplot_start):
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='massdiff_ppm'), fig, subplot_max_x, subplot_start, xlabel='precursor mass difference, ppm')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='RT diff'), fig, subplot_max_x, subplot_start, xlabel='RT difference, min')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='num_missed_cleavages'), fig, subplot_max_x, subplot_start, xlabel='missed cleavages')
    subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='assumed_charge'), fig, subplot_max_x, subplot_start, xlabel='precursor charge')
    subplot_start += 1
    if 'fragmentMT' in df.columns:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='fragmentMT'), fig, subplot_max_x, subplot_start, xlabel='median fragment error, Da')
        subplot_start += 1
    if 'ISOWIDTHDIFF' in df.columns:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='ISOWIDTHDIFF'), fig, subplot_max_x, subplot_start, xlabel='isolation mass error, Da')
        subplot_start += 1
    if 'MS1Intensity' in df.columns:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='MS1Intensity'), fig, subplot_max_x, subplot_start, xlabel='MS1 Intensity')
        subplot_start += 1

    if len(set(df['massdiff_int'])) > 1:
        plot_hist_descriptor(get_descriptor_array(df, df_f, dname='massdiff_int'), fig, subplot_max_x, subplot_start, xlabel='isotope mass difference, Da')
        subplot_start += 1
    for df_col in df.columns:
        if df_col.startswith('mass shift'):
            plot_hist_descriptor(get_descriptor_array(df, df_f, dname=df_col), fig, subplot_max_x, subplot_start, xlabel=df_col)
            subplot_start += 1
    plot_hist_descriptor(get_descriptor_array(df, df_f, dname='log_score'), fig, subplot_max_x, subplot_start, xlabel='LOG10(ML score)')
    subplot_start += 1
    separate_figures = isinstance(fig, str)
    if not separate_figures:
        plot_legend(fig, subplot_max_x, subplot_start)
    subplot_start += 1

def get_bins(inarrays, bin_size_one=False):
    tmp = np.concatenate(inarrays)
    minv = tmp.min()
    maxv = tmp.max()
    if bin_size_one:
        return np.arange(minv, maxv+1, 1)
    else:
        return np.linspace(minv, maxv+1, num=100)

def get_bins_for_descriptors(inarrays):
    tmp = np.concatenate(inarrays)
    minv = tmp.min()
    maxv = tmp.max()
    if len(set(tmp)) <= 15:
        return np.arange(minv, maxv + 2, 1.0), 1.0
    binsize = False
    for inar in inarrays:
        binsize_tmp = get_fdbinsize(inar)
        if not binsize or binsize > binsize_tmp:
            binsize = binsize_tmp
    # binsize = get_fdbinsize(tmp)
    if binsize < float(maxv - minv) / 300:
        binsize = float(maxv - minv) / 300
    lbin_s = scoreatpercentile(tmp, 1.0)
    lbin = minv
    if lbin_s and abs((lbin - lbin_s) / lbin_s) > 1.0:
        lbin = lbin_s * 1.05
    rbin_s = scoreatpercentile(tmp, 99.0)
    rbin = maxv
    if rbin_s and abs((rbin - rbin_s) / rbin_s) > 1.0:
        rbin = rbin_s * 1.05
    rbin += 1.5 * binsize
    return np.arange(lbin, rbin + binsize, binsize), binsize

def get_fdbinsize(data_list):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    data_list:  1D Data set
    Returns:
    optimal_bin_size:  F-D bin size
    """
    data_list = np.sort(data_list)
    upperquartile = scoreatpercentile(data_list, 75)
    lowerquartile = scoreatpercentile(data_list, 25)
    iqr = upperquartile - lowerquartile
    optimal_bin_size = 2. * iqr / len(data_list) ** (1. / 3.)
    return optimal_bin_size

def normalize_fname(s):
    return re.sub(r'[<>:\|/?*]', '', s)

def outpath(outfolder, s, ext='.png'):
    return os.path.join(outfolder, normalize_fname(s) + ext)

def plot_outfigures(df, df_f, df_peptides, df_peptides_f, outfolder, outbasename, df_proteins, df_proteins_f, separate_figures=False):
    if not separate_figures:
        fig = plt.figure(figsize=(16, 12))
        dpi = fig.get_dpi()
        fig.set_size_inches(3000.0/dpi, 3000.0/dpi)
    else:
        outfolder = os.path.join(outfolder, outbasename + '_figures')
        os.makedirs(outfolder, exist_ok=True)
        fig = outfolder
    subplot_max_x = calc_max_x_value(df, df_proteins)
    descriptor_start_index = 7
    plot_basic_figures(df, df_f, fig, subplot_max_x, 1, 'PSMs')
    plot_basic_figures(df_peptides, df_peptides_f, fig, subplot_max_x, 4, 'peptides')
    if 'LOG10_NSAF' in df_proteins.columns:
        plot_protein_figures(df_proteins, df_proteins_f, fig, subplot_max_x, 7)
        plot_aa_stats(df_f, df_proteins_f, fig, subplot_max_x, 9)
        descriptor_start_index += 3
    plot_descriptors_figures(df, df_f, fig, subplot_max_x, descriptor_start_index)
    plt.grid(color='#EEEEEE')
    plt.tight_layout()
    if not separate_figures:
        plt.savefig(os.path.join(outfolder, outbasename) + '.pdf')
