import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome import refgenome
from sklearn import preprocessing
from statannot import add_stat_annotation
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('-ic', '--input_cn', type=str, nargs='+', help='list of pseudobulk CN profiles for each sample')
    p.add_argument('-ir', '--input_rt', type=str, nargs='+', help='list of pseudobulk RT profiles for each sample')
    p.add_argument('-sc', '--sample_corrs', type=str, help='plot of CN and RT correlations between samples')
    p.add_argument('-cc', '--clone_corrs', type=str, help='plot of CN and RT correlations between clones')

    return p.parse_args()

def read_rt_data(argv):
    # read in the pseudobulk RT profiles for each sample
    rt = pd.DataFrame()
    for path in argv.input_rt:
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

    return rt


def read_cn_data(argv):
    # read in the pseudobulk CN profiles for each sample
    cn = pd.DataFrame()

    for file in argv.input_cn:
        temp_cn = pd.read_csv(file, sep='\t')

        # find the dataset name given the path
        d = file.split('/')[2]

        # rename the 'sample_{dataset}' column to '{dataset}_sample_cn'
        if 'sample_{}'.format(d) in temp_cn.columns:
            temp_cn.rename(columns={'sample_{}'.format(d): '{}_sample_cn'.format(d)}, inplace=True)
        # sometimes the sample column is named 'dataset_{dataset}'
        elif 'dataset_{}'.format(d) in temp_cn.columns:
            temp_cn.rename(columns={'dataset_{}'.format(d): '{}_sample_cn'.format(d)}, inplace=True)

        # add the dataset name as prefix to all columns with 'clone' in the name
        cols = [c for c in temp_cn.columns if c.startswith('clone')]
        for c in cols:
            temp_cn.rename(columns={c: '{}_{}_cn'.format(d, c)}, inplace=True)
        
        # merge the temp_cn dataframe with the cn dataframe if cn is not empty
        if cn.empty:
            cn = temp_cn
        else:
            cn = pd.merge(cn, temp_cn)

    # set chr column to category
    cn.chr = cn.chr.astype('str')
    cn.chr = cn.chr.astype('category')

    cn['end'] = cn['start'] + 500000 - 1

    return cn


def plot_sample_corrs(cn, rt, sample_cn_cols, sample_rt_cols, argv):
    ''' Create a figure that shows the correlation between the sample CN and RT profiles. '''
    fig, ax = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)

    # correlation between copy number profiles
    cn_corrs = cn[sample_cn_cols].corr()
    mask = np.zeros_like(cn_corrs, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(cn_corrs, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='.2f', ax=ax[0])
    ax[0].set_title('Correlation in sample CN')

    # correlation between RT profiles
    rt_corrs = rt[sample_rt_cols].corr()
    mask = np.zeros_like(rt_corrs, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(rt_corrs, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='.2f', ax=ax[1])
    ax[1].set_title('Correlation in sample RT')

    # only include the sample_id prefix in the xticklabels and yticklabels
    for a in ax:
        a.set_xticklabels([str(c.get_text()).split('_')[0] for c in a.get_xticklabels()])
        a.set_yticklabels([str(c.get_text()).split('_')[0] for c in a.get_yticklabels()])

    fig.savefig(argv.sample_corrs, dpi=300, bbox_inches='tight')


def plot_clone_corrs(cn, rt, clone_cn_cols, clone_rt_cols, argv):
    ''' Create figure that shows the correlation between the clone CN and RT profiles across all samples. '''
    fig, ax = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)

    # correlation between copy number profiles
    cn_corrs = cn[clone_cn_cols].corr()
    mask = np.zeros_like(cn_corrs, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(cn_corrs, square=False, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, ax=ax[0])
    ax[0].set_title('Correlation in clone CN')

    # correlation between RT profiles
    rt_corrs = rt[clone_rt_cols].corr()
    mask = np.zeros_like(rt_corrs, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(rt_corrs, square=False, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, ax=ax[1])
    ax[1].set_title('Correlation in clone RT')

    # only include the sample_id prefix and clone_id letter in the xticklabels and yticklabels
    ax[0].set_xticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2]) for c in ax[0].get_xticklabels()])
    ax[0].set_yticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2]) for c in ax[0].get_yticklabels()])

    # naming convention for the clone RT profiles is different than the clone CN profiles
    ax[1].set_xticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2].replace('clone', '')) for c in ax[1].get_xticklabels()])
    ax[1].set_yticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2].replace('clone', '')) for c in ax[1].get_yticklabels()])
 
    fig.savefig(argv.clone_corrs, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    # read in the pseudobulk RT profiles for all samples
    rt = read_rt_data(argv)

    # read in the pseudobulk CN profiles for all samples
    cn = read_cn_data(argv)

    # create a list of the sample cn profiles
    sample_cn_cols = [c for c in cn.columns if c.endswith('sample_cn')]
    # create a list of the sample rt profiles
    sample_rt_cols = [c for c in rt.columns if c.endswith('model_rep_state') and 'clone' not in c]
    # create a list of the clone cn profiles
    clone_cn_cols = [c for c in cn.columns if c.endswith('cn') and 'sample' not in c]
    # create a list of the clone rt profiles
    clone_rt_cols = [c for c in rt.columns if c.endswith('model_rep_state') and 'clone' in c]

    # plot the correlation between the sample CN and RT profiles
    plot_sample_corrs(cn, rt, sample_cn_cols, sample_rt_cols, argv)

    # plot the correlation between the clone CN and RT profiles
    plot_clone_corrs(cn, rt, clone_cn_cols, clone_rt_cols, argv)


if __name__ == '__main__':
    main()
