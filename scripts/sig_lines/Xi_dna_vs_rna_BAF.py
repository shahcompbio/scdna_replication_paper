import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_htert_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('--datasets', type=str, nargs='+', help='name of each dataset for extracting SIGNALS results')
    p.add_argument('--table', help='table of per-chromosome DNA vs RNA BAFs across all datasets')
    p.add_argument('--plot', help='Plot of DNA vs RNA BAFs across all chromosomes and all datasets')

    return p.parse_args()


# load scRNA signals results for a test sample
def load_scrna(dataset):
    path = '/juno/work/shah/users/william1/projects/ascn_v2/results/scrna/counthaps/{d}/{d}-allele-counts-phased.csv.gz'.format(d=dataset)
    df = pd.read_csv(path, dtype={'chr': str})
    df['dataset'] = dataset
    return df

# load the scDNA signals results
def load_scdna(dataset):
    path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/schnapps-results/persample/{d}_hscn.csv.gz'.format(d=dataset)
    df = pd.read_csv(path, dtype={'chr': str})
    return df


def compute_chrom_baf_df(df, modality='rna'):
    '''
    Given a dataframe of SIGNALS results, compute the mean and std BAF for each chromosome across all cells.
    '''
    baf_col = '{}_baf'.format(modality)
    baf_mean_col = '{}_baf_mean'.format(modality)
    baf_std_col = '{}_baf_std'.format(modality)
    # count the total number of A and B allele counts per cell per chromosome
    df = df[['cell_id', 'chr', 'patient', 'alleleA', 'alleleB']].groupby(['patient', 'cell_id', 'chr']).sum().reset_index()
    # compute the RNA BAF within each cell & chromosome
    df[baf_col] = df['alleleB'] / (df['alleleA'] + df['alleleB'])
    # compute the mean RNA BAF for each chromosome across all cells
    df[baf_mean_col] = df.groupby(['chr'])[baf_col].transform('mean')
    # compute the standard deviation of RNA BAF for each chromosome across all cells
    df[baf_std_col] = df.groupby(['chr'])[baf_col].transform('std')
    # subset to just the columns of interest
    df = df[['patient', 'chr', baf_mean_col, baf_std_col]].drop_duplicates()
    # add a column to indicate whether the bin is in chrX or an autosome
    df['chr_type'] = np.where(df['chr'] == 'X', 'chrX', 'autosome')
    return df


def main():
    argv = get_args()

    df = []
    for d in argv.datasets:
        temp_rna = compute_chrom_baf_df(load_scrna(d), modality='rna')
        temp_dna = compute_chrom_baf_df(load_scdna(d), modality='dna')
        temp_df = pd.merge(temp_rna, temp_dna)
        df.append(temp_df)
    df = pd.concat(df, ignore_index=True)

    # save the table
    df.to_csv(argv.table, index=False)

    # plot the table
    # points are colored only for chrX and each sample is a different color
    # autosomes are grey
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    sns.scatterplot(data=df.query("chr_type=='autosome'"), x='dna_baf_mean', y='rna_baf_mean', color='grey', alpha=0.5, ax=ax)
    sns.scatterplot(data=df.query("chr_type=='chrX'"), x='dna_baf_mean', y='rna_baf_mean', hue='patient', palette=get_htert_cmap(), ax=ax)
    ax.plot([0, 1], [0, 1], linestyle='--', color='black')
    ax.set_xlabel('DNA BAF')
    ax.set_ylabel('RNA BAF')
    ax.set_title('hTERT cell lines\nMean BAF per chromosome per sample')
    ax.legend(title='sample chrX', loc='lower right')
    fig.savefig(argv.plot, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
