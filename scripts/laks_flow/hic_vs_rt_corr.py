import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('hic', help='HiC compartments for select cell lines')
    p.add_argument('rt', help='pseudobulk RT profiles from T47D and GM18507 cell lines')
    p.add_argument('output', help='plot of RT vs. HiC compartment pairwise correlations')

    return p.parse_args()


def main():
    argv = get_args()

    hic = pd.read_csv(argv.hic)
    rt = pd.read_csv(argv.rt, sep='\t')

    print('rt.head():', rt.head(), sep='\n')

    # drop columns that contain the 'merged' or 'diff' substrings
    rt = rt.loc[:, ~rt.columns.str.contains('merged|diff')]
    # remove the '_split' substring
    rt.columns = rt.columns.str.replace('_split', '')

    print('hic.head():', hic.head(), sep='\n')
    print('rt.head():', rt.head(), sep='\n')

    # create one dataframe that has hic and rt data merged
    hic_rt = pd.merge(hic, rt)
    print('hic_rt.head():', hic_rt.head(), sep='\n')

    # subset rt columns to just the representative cell lines
    rt_cols = ['rt_GM18507', 'rt_T47D']

    # subset hic columns to just the representative cell lines
    hic_cols = [
        'GM11168_HiC_compartment_score', 'GM12878_HiC_compartment_score', 
        'GM13976_HiC_compartment_score', 'GM13977_HiC_compartment_score', 
        'GM18951_HiC_compartment_score', 'T47D_HiC_compartment_score_hg19',
    ]

    corr_mat = pd.DataFrame(index=hic_cols, columns=rt_cols)
    for hic_col in hic_cols:
        for rt_col in rt_cols:
            # drop the NaN values for these two columns before taking the correlation
            temp = hic_rt[[hic_col, rt_col]].dropna()
            # find the correlation between the two columns
            corr, pval = spearmanr(temp[hic_col], temp[rt_col])
            # add the corr to corr_mat dataframe
            corr_mat.loc[hic_col, rt_col] = float(corr)
    
    print('corr_mat:', corr_mat, sep='\n')
    
    # rename the columns to drop the 'rt_' substring
    corr_mat.columns = corr_mat.columns.str.replace('rt_', '')

    # rename the rows to only use te cell line prefix before the '_HiC' substring
    corr_mat.index = corr_mat.index.str.replace('_HiC_compartment_score', '')
    corr_mat.index = corr_mat.index.str.replace('_hg19', '')
    corr_mat.index = corr_mat.index.str.replace('MammEp', 'Mammary Epithelial')

    print('corr_mat:', corr_mat, sep='\n')

    # convert corr_mat to dtype float
    corr_mat = corr_mat.astype(float)

    print('corr_mat:', corr_mat, sep='\n')

    # plot as a seaborn heatmap
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    sns.heatmap(corr_mat, square=False, linewidths=.5, cbar_kws={"shrink": .5}, annot=True, fmt='.2f', ax=ax)
    ax.set_xlabel('PERT pseudobulk RT')
    ax.set_ylabel('Hi-C from ENCODE')
    ax.set_title('RT vs Hi-C compartment\nSpearman correlation')
    # center the yticklabels at the middle of the tick
    ax.set_yticklabels(ax.get_yticklabels(), va='center')

    fig.savefig(argv.output, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
