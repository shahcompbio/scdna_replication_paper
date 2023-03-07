import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome import refgenome
from sklearn import preprocessing
from statannot import add_stat_annotation
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from common.plot_utils import plot_cell_cn_profile2


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='list of pseudobulk RT profiles for each sample')
    p.add_argument('-c', '--counts', type=str, help='cell cycle clone counts for across all samples')
    p.add_argument('-srp', '--sample_rt_profiles', type=str, help='plot of the RT profiles for each sample')
    p.add_argument('-srd', '--sample_rt_diffs', type=str, help='violin plots of difference in RT between autosomes and chrX for samples with and without chrX CNAs')
    p.add_argument('-crp', '--clone_rt_profiles', type=str, help='plot of the RT profiles for each clone, one subpanel per sample')
    p.add_argument('-crd4', '--clone_rt_diffs_SA1054', type=str, help='violin plots of difference in RT between autosomes and chrX for clones with and without chrX CNAs for SA1054')
    p.add_argument('-crd5', '--clone_rt_diffs_SA1055', type=str, help='violin plots of difference in RT between autosomes and chrX for clones with and without chrX CNAs for SA1055')


    return p.parse_args()


def violins_with_pvals(df, x, y, hue, ax, box_pairs, order=None, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, order=order)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test, order=order,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_rt_diff_vs_chrX(df, ax, x='chrX', y='RT_diff', test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' 
    Plot the difference in pseudobulk RT values between SA039 and OV2295
    where the data is split by chrX vs autosomes
    '''
    x = x
    y = y
    hue = None
    box_pairs = [
        ('autosomes', 'chrX'),
    ]
    order = ['autosomes', 'chrX']
    violins_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose)
    return ax


def plot_sample_rt_profiles(rt, rt_coi, counts, argv):
    fig, ax = plt.subplots(2, 1, figsize=(16,8), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_coi):
        dataset_id = col.split('_')[0]

        # find the number of S-phase cells in this dataset according to counts
        n = int(np.sum(counts.query('dataset=="{}"'.format(dataset_id))['num_cells_s'].values))
        label = '{} (n={})'.format(dataset_id, n)
        
        # plot pseudobulk rt profile for this dataset
        plot_cell_cn_profile2(
            ax[0], rt, col, color='C{}'.format(i), 
            max_cn=None, scale_data=False, lines=True, label=label
        )
        # zoom in on chrX
        plot_cell_cn_profile2(
            ax[1], rt, col, color='C{}'.format(i), chromosome='X',
            max_cn=None, scale_data=False, lines=True, label=label
        )

    ax[0].set_title('Sample Pseudobulk RT')
    ax[1].set_title('Sample Pseudobulk RT')
    ax[0].set_ylabel('RT profile\n<--late | early-->')
    ax[1].set_ylabel('RT profile\n<--late | early-->')
    ax[0].legend(title='Sample ID', loc='upper left')
    ax[1].legend(title='Sample ID', loc='upper left')

    fig.savefig(argv.sample_rt_profiles, dpi=300, bbox_inches='tight')


def plot_sample_rt_diffs(rt_diff, rt_diff_coi, argv):
    fig, ax = plt.subplots(2, 2, figsize=(6,8), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_diff_coi):
        d = col.split('_')[0]
        plot_rt_diff_vs_chrX(rt_diff, ax[i], y=col)
        ax[i].set_xlabel(None)
        ax[i].set_ylabel('Difference in pseudobulk RT\n<--ref earlier | {} earlier -->'.format(d))
        ax[i].set_title('chrX RT in {}'.format(d))

    fig.savefig(argv.sample_rt_diffs, dpi=300, bbox_inches='tight')



def plot_clone_rt_profiles(rt, ax0, ax1, dataset_id, counts):
    # get the clone RT columns for this dataset
    clone_rt_cols = [c for c in rt.columns if c.startswith(dataset_id) and 'clone' in c]

    # plot each clone RT profile as a different color
    # only plot the X chromosome for now
    for i, col in enumerate(clone_rt_cols):
        clone_id = col.split('_')[2].split('clone')[1]
        # find the number of S-phase cells in this clone according to dataset_counts
        n = int(counts.query("dataset=='{}'".format(dataset_id)).query("clone_id=='{}'".format(clone_id))['num_cells_s'].values[0])
        # plot the X chromosome
        plot_cell_cn_profile2(
            ax1, rt, col, color='C{}'.format(i), chromosome='X',
            max_cn=None, scale_data=False, lines=True, label='{} (n={})'.format(clone_id, n)
        )
        # plot the whole genome
        plot_cell_cn_profile2(
            ax0, rt, col, color='C{}'.format(i),
            max_cn=None, scale_data=False, lines=True, label='{} (n={})'.format(clone_id, n)
        )
    
    # format the plot
    for ax in [ax0, ax1]:
        ax.set_title('{} clone RT profiles'.format(dataset_id))
        ax.set_ylabel('RT profile\n<--late | early-->')
        ax.legend(title='Clone ID', loc='upper left')


def plot_clone_rt_profiles_wrapper(rt, counts, argv):
    # plot clone RT profiles for each dataset
    # find the number of datasets
    datasets = set([c.split('_')[0] for c in rt.columns if (c.startswith('SA') or c.startswith('OV')) and 'clone' not in c])
    n_datasets = len(datasets)

    # create a grid of subplots where each row is a dataset and the height is propotional to the number of datasets
    fig, ax = plt.subplots(n_datasets, 2, figsize=(32, 4*n_datasets), tight_layout=True)

    # loop through each dataset and plot the clone RT profiles
    for i, dataset_id in enumerate(datasets):
        plot_clone_rt_profiles(rt, ax[i, 0], ax[i, 1], dataset_id, counts)
    
    # save the figure
    fig.savefig(argv.clone_rt_profiles, dpi=300, bbox_inches='tight')


def plot_clone_rt_diffs_SA1054(rt, counts, argv):
    ''' compute RT difference between the SA1054 clones based on allelic state '''
    # list of clones that are allelically balanced
    ref_rt_clones_SA1054 = [
        'SA1054_pseudobulk_cloneB_model_rep_state',
        'SA1054_pseudobulk_cloneC_model_rep_state',
    ]

    # use the number of cells in each clone to compute a corresponding weight vector for the reference
    ref_rt_weights_SA1054 = [
        counts.query("dataset=='SA1054'").query("clone_id=='B'")['num_cells_s'].values[0],
        counts.query("dataset=='SA1054'").query("clone_id=='C'")['num_cells_s'].values[0],
    ]

    # combine clones D & E into one pseudobulk RT profile since they have the same allelic state
    # find the number of cells in clones D & E
    SA1054_clones_DE = [
        'SA1054_pseudobulk_cloneD_model_rep_state',
        'SA1054_pseudobulk_cloneE_model_rep_state',
    ]
    weights_SA1054_clones_DE = [
        counts.query("dataset=='SA1054'").query("clone_id=='D'")['num_cells_s'].values[0],
        counts.query("dataset=='SA1054'").query("clone_id=='E'")['num_cells_s'].values[0],
    ]

    # take a weighted mean of the clone D & E pseudobulk RT profiles
    rt['SA1054_pseudobulk_cloneDE_model_rep_state'] = np.average(rt[SA1054_clones_DE], weights=weights_SA1054_clones_DE, axis=1)

    # list of clones that are allelically imbalanced
    rt_coi_clones_SA1054 = [
        'SA1054_pseudobulk_cloneA_model_rep_state',
        'SA1054_pseudobulk_cloneDE_model_rep_state',
    ]

    # subset the RT dataframe to only include the clones of interest or the ref clones
    rt_coi_SA1054 = ref_rt_clones_SA1054 + rt_coi_clones_SA1054 + ['chr', 'start', 'end']
    rt_diff_SA1054 = rt[rt_coi_SA1054]

    # denote which loci are part of chrX
    rt_diff_SA1054['chrX'] = rt_diff_SA1054['chr'].apply(lambda x: 'chrX' if x=='X' else 'autosomes')

    # use the ref clones to compute the reference RT
    rt_diff_SA1054['ref_rt'] = np.average(rt_diff_SA1054[ref_rt_clones_SA1054], weights=ref_rt_weights_SA1054, axis=1)

    # for all the clones in rt_coi_clones_SA1054, compute the RT difference between the clone and the ref_rt column
    rt_diff_columns_SA1054 = []
    for col in rt_coi_clones_SA1054:
        clone_id = col.split('_')[2].split('clone')[1]
        temp_diff_col = '{}_rt_diff'.format(clone_id)
        rt_diff_SA1054[temp_diff_col] = rt_diff_SA1054[col] - rt_diff_SA1054['ref_rt']
        rt_diff_columns_SA1054.append(temp_diff_col)

    # for every rt_diff_column, plot the RT difference vs chrX
    fig, ax = plt.subplots(1, 2, figsize=(6,4), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_diff_columns_SA1054):
        clone_id = col.split('_')[0]
        plot_rt_diff_vs_chrX(rt_diff_SA1054, ax[i], y=col)
        ax[i].set_xlabel(None)
        ax[i].set_ylabel('Difference in pseudobulk RT\n<--ref earlier | {} earlier -->'.format(clone_id))
        ax[i].set_title('chrX RT in SA1054 clone {}'.format(clone_id))

    fig.savefig(argv.clone_rt_diffs_SA1054, dpi=300, bbox_inches='tight')


def plot_clone_rt_diffs_SA1055(rt, counts, argv):
    ''' compute RT difference between the SA1055 clones based on allelic state '''
    # list of clones that are allelically balanced
    ref_rt_clones_SA1055 = [
        'SA1055_pseudobulk_cloneE_model_rep_state',
        'SA1055_pseudobulk_cloneF_model_rep_state',
        'SA1055_pseudobulk_cloneG_model_rep_state',
        'SA1055_pseudobulk_cloneH_model_rep_state',
        'SA1055_pseudobulk_cloneI_model_rep_state',
    ]

    # use the number of cells in each clone to compute a corresponding weight vector for the reference
    ref_rt_weights_SA1055 = [
        counts.query("dataset=='SA1055'").query("clone_id=='E'")['num_cells_s'].values[0],
        counts.query("dataset=='SA1055'").query("clone_id=='F'")['num_cells_s'].values[0],
        counts.query("dataset=='SA1055'").query("clone_id=='G'")['num_cells_s'].values[0],
        counts.query("dataset=='SA1055'").query("clone_id=='H'")['num_cells_s'].values[0],
        counts.query("dataset=='SA1055'").query("clone_id=='I'")['num_cells_s'].values[0],
    ]

    # combine clones A & B into one pseudobulk RT profile since they have the same allelic state
    # find the number of cells in clones D & E
    SA1055_clones_AB = [
        'SA1055_pseudobulk_cloneA_model_rep_state',
        'SA1055_pseudobulk_cloneB_model_rep_state',
    ]
    weights_SA1055_clones_AB = [
        counts.query("dataset=='SA1055'").query("clone_id=='A'")['num_cells_s'].values[0],
        counts.query("dataset=='SA1055'").query("clone_id=='B'")['num_cells_s'].values[0],
    ]

    # take a weighted mean of the clone A & B pseudobulk RT profiles
    rt['SA1055_pseudobulk_cloneAB_model_rep_state'] = np.average(rt[SA1055_clones_AB], weights=weights_SA1055_clones_AB, axis=1)

    # list of clones that are allelically imbalanced
    rt_coi_clones_SA1055 = [
        'SA1055_pseudobulk_cloneAB_model_rep_state',
        'SA1055_pseudobulk_cloneC_model_rep_state',
        'SA1055_pseudobulk_cloneD_model_rep_state',
    ]

    # subset the RT dataframe to only include the clones of interest or the ref clones
    rt_coi_SA1055 = ref_rt_clones_SA1055 + rt_coi_clones_SA1055 + ['chr', 'start', 'end']
    rt_diff_SA1055 = rt[rt_coi_SA1055]

    # denote which loci are part of chrX
    rt_diff_SA1055['chrX'] = rt_diff_SA1055['chr'].apply(lambda x: 'chrX' if x=='X' else 'autosomes')

    # use the ref clones to compute the reference RT
    rt_diff_SA1055['ref_rt'] = np.average(rt_diff_SA1055[ref_rt_clones_SA1055], weights=ref_rt_weights_SA1055, axis=1)

    # for all the clones in rt_coi_clones_SA1055, compute the RT difference between the clone and the ref_rt column
    rt_diff_columns_SA1055 = []
    for col in rt_coi_clones_SA1055:
        clone_id = col.split('_')[2].split('clone')[1]
        temp_diff_col = '{}_rt_diff'.format(clone_id)
        rt_diff_SA1055[temp_diff_col] = rt_diff_SA1055[col] - rt_diff_SA1055['ref_rt']
        rt_diff_columns_SA1055.append(temp_diff_col)

    # for every rt_diff_column, plot the RT difference vs chrX
    fig, ax = plt.subplots(1, 3, figsize=(9,4), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_diff_columns_SA1055):
        clone_id = col.split('_')[0]
        plot_rt_diff_vs_chrX(rt_diff_SA1055, ax[i], y=col)
        ax[i].set_xlabel(None)
        ax[i].set_ylabel('Difference in pseudobulk RT\n<--ref earlier | {} earlier -->'.format(clone_id))
        ax[i].set_title('chrX RT in SA1055 clone {}'.format(clone_id))

    fig.savefig(argv.clone_rt_diffs_SA1055, dpi=300, bbox_inches='tight')


def main():
    argv = get_args()

    # load the cell counts in each clone
    counts = pd.read_csv(argv.counts, sep='\t')

    # read in the pseudobulk RT profiles for each sample
    rt = pd.DataFrame()
    for path in argv.input:
        temp_rt = pd.read_csv(path, sep='\t')
        cols = [c for c in temp_rt.columns if c.startswith('pseudobulk') and c.endswith('model_rep_state')]
        cols.append('chr')
        cols.append('start')
        temp_rt = temp_rt[cols]
        
        # add dataset name as prefix to RT columns
        d = path.split('/')[2]
        for c in temp_rt.columns:
            if c.startswith('pseudobulk') and c.endswith('model_rep_state'):
                temp_rt.rename(columns={c: '{}_{}'.format(d, c)}, inplace=True)
        
        if rt.empty:
            rt = temp_rt
        else:
            rt = pd.merge(rt, temp_rt)

    # set chr column to category
    rt.chr = rt.chr.astype('str')
    rt.chr = rt.chr.astype('category')

    # add end position as it's necessary for plotting functions
    rt['end'] = rt['start'] + 500000 - 1

    # columns of interest for sample RT profiles
    rt_coi = [c for c in rt.columns if (c.endswith('_pseudobulk_model_rep_state'))]

    # plot the sample RT pseudobulk profiles
    plot_sample_rt_profiles(rt, rt_coi, counts, argv)

    # samples with no CNAs on chrX should be used as a reference
    ref_rt_samples = [
        'SA039_pseudobulk_model_rep_state',
        'SA906a_pseudobulk_model_rep_state',
        'SA906b_pseudobulk_model_rep_state',
        'SA1292_pseudobulk_model_rep_state',
        'SA1188_pseudobulk_model_rep_state',
    ]
    ref_rt_weights = [
        np.sum(counts.query('dataset=="SA039"')['num_cells_s'].values),
        np.sum(counts.query('dataset=="SA906a"')['num_cells_s'].values),
        np.sum(counts.query('dataset=="SA906b"')['num_cells_s'].values),
        np.sum(counts.query('dataset=="SA1292"')['num_cells_s'].values),
        np.sum(counts.query('dataset=="SA1188"')['num_cells_s'].values),
    ]

    # samples with CNAs on chrX should be used as a test
    test_rt_samples = [
        'SA1056_pseudobulk_model_rep_state',
        'SA1054_pseudobulk_model_rep_state',
        'SA1055_pseudobulk_model_rep_state',
        'OV2295_pseudobulk_model_rep_state'
    ]

    # create a new dataframe with the difference in RT values between the reference samples
    rt_diff = rt[[c for c in rt.columns if (c in rt_coi or c in ['chr', 'start', 'end'])]]

    # use the average of all the reference samples as the reference RT profile
    rt_diff['ref_rt'] = np.average(rt[ref_rt_samples], weights=ref_rt_weights, axis=1)

    # denote which loci are on chrX
    rt_diff['chrX'] = rt_diff['chr'].apply(lambda x: 'chrX' if x=='X' else 'autosomes')

    # calculate the difference in RT values between the reference and test samples
    for col in test_rt_samples:
        dataset_id = col.split('_')[0]
        rt_diff['{}_rt_diff'.format(dataset_id)] = rt_diff[col] - rt_diff['ref_rt']
    
    # denote which columns are of interest for plotting
    rt_diff_coi = [
        'SA1056_rt_diff',
        'SA1054_rt_diff',
        'SA1055_rt_diff',
        'OV2295_rt_diff'
    ]

    # plot the difference in RT values between the reference and test samples
    plot_sample_rt_diffs(rt_diff, rt_diff_coi, argv)

    # plot the clone RT profiles for chrX
    plot_clone_rt_profiles_wrapper(rt, counts, argv)

    # plot clone-specific RT differences on chrX for SA1054 and SA1055
    plot_clone_rt_diffs_SA1054(rt, counts, argv)
    plot_clone_rt_diffs_SA1055(rt, counts, argv)


if __name__ == '__main__':
    main()
