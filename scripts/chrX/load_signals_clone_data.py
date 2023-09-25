import pandas as pd
import numpy as np
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, help='chrX RT delays for clone profiles')
    p.add_argument('-o', '--output', type=str, help='SIGNALS DNA and RNA BAF for each chromosome arm')

    return p.parse_args()


def add_p_q_arm(df):
    centromere_dict = {
        '1': 121535434,
        '2': 92326171,
        '3': 90504854,
        '4': 49660117,
        '5': 46405641,
        '6': 58830166,
        '7': 58054331,
        '8': 43838887,
        '9': 47367679,
        '10': 39254935,
        '11': 51644205,
        '12': 34856694,
        '13': 16000000,
        '14': 16000000,
        '15': 17000000,
        '16': 35335801,
        '17': 22263006,
        '18': 15460898,
        '19': 24681782,
        '20': 26369569,
        '21': 11288129,
        '22': 13000000,
        'X': 58632012,
        'Y': 57217415
    }

    for chrom, chr_df in df.groupby('chr'):
        if chrom in centromere_dict:
            centromere = centromere_dict[chrom]
            df.loc[((df['chr']==chrom) & (df['start'] < centromere)), 'p_q_arm'] = 'p'
            df.loc[((df['chr']==chrom) & (df['start'] >= centromere)), 'p_q_arm'] = 'q'
    df['chr_arm'] = df['chr'] + df['p_q_arm']
    
    return df


# load scRNA signals results for a test sample
def load_scrna(dataset):
    path = '/juno/work/shah/users/william1/projects/ascn_v2/results/scrna/counthaps/{d}/{d}-allele-counts-phased.csv.gz'.format(d=dataset)
    df = pd.read_csv(path, dtype={'chr': str})
    df['dataset'] = dataset
    df = add_p_q_arm(df)
    return df

# load the scDNA signals results with clone_id labels merged in
def load_scdna(dataset):
    path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/schnapps-results/persample/{d}_hscn.csv.gz'.format(d=dataset)
    df = pd.read_csv(path, dtype={'chr': str})
    # load reclustered data for fitness datasets
    if dataset in ['SA535', 'SA1035', 'SA609']:
        clones = pd.read_csv('/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/fitness/{d}/cn_data_clones.tsv'.format(d=dataset), sep='\t')
        clones = clones[['cell_id', 'clone_id']].drop_duplicates()
    else:
        # otherwise load ths SIGNALS clustering
        clones = pd.read_csv('/juno/work/shah/users/weinera2/projects/scdna_replication_paper/data/signatures/clone_trees/{d}_clones.tsv'.format(d=dataset), sep='\t')
    df = pd.merge(df, clones)
    df = add_p_q_arm(df)
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
    # add a dummy column for chr_arm ('f' for full chromosome)
    df['chr_arm'] = df['chr'] + 'f'
    return df


def compute_arm_baf_df(df, modality='rna'):
    '''
    Given a dataframe of SIGNALS results, compute the mean and std BAF for each chromosome arm across all cells.
    '''
    baf_col = '{}_baf'.format(modality)
    baf_mean_col = '{}_baf_mean'.format(modality)
    baf_std_col = '{}_baf_std'.format(modality)
    # count the total number of A and B allele counts per cell per chromosome arm
    df = df[['cell_id', 'chr', 'chr_arm', 'patient', 'alleleA', 'alleleB']].groupby(['patient', 'cell_id', 'chr', 'chr_arm']).sum().reset_index()
    # compute the RNA BAF within each cell & chromosome
    df[baf_col] = df['alleleB'] / (df['alleleA'] + df['alleleB'])
    # compute the mean RNA BAF for each chromosome across all cells
    df[baf_mean_col] = df.groupby(['chr_arm'])[baf_col].transform('mean')
    # compute the standard deviation of RNA BAF for each chromosome across all cells
    df[baf_std_col] = df.groupby(['chr_arm'])[baf_col].transform('std')
    # subset to just the columns of interest
    df = df[['patient', 'chr', 'chr_arm', baf_mean_col, baf_std_col]].drop_duplicates()
    # add a column to indicate whether the bin is in chrX or an autosome
    df['chr_type'] = np.where(df['chr'] == 'X', 'chrX', 'autosome')
    return df


def main():
    argv = get_args()

    # load the input mean rt file that contains sample RT delays
    mean_rt = pd.read_csv(argv.input)

    datasets = mean_rt['dataset'].unique().tolist()
    df = []
    for d in datasets:
        # load in the scDNA SIGNALS data with clone_id labels merged in
        temp_scdna = load_scdna(d)
        for c in temp_scdna.clone_id.unique():
            # subset to just the clone of interest
            temp_scdna_clone = temp_scdna.query('clone_id==@c')
            # compute the mean and std BAF for each chromosome across all cells
            temp_scdna_clone_baf_chr = compute_chrom_baf_df(temp_scdna_clone, modality='dna')
            # add a column to indicate the clone_id
            temp_scdna_clone_baf_chr['clone_id'] = c
            # add a column to indicate the dataset
            temp_scdna_clone_baf_chr['dataset'] = d
            # compute the mean and std for each chromsome arm across all cells
            temp_scdna_clone_baf_arm = compute_arm_baf_df(temp_scdna_clone, modality='dna')
            # add a column to indicate the clone_id
            temp_scdna_clone_baf_arm['clone_id'] = c
            # add a column to indicate the dataset
            temp_scdna_clone_baf_arm['dataset'] = d
            # concatenate the two dataframes
            temp_scdna_clone_baf = pd.concat([temp_scdna_clone_baf_chr, temp_scdna_clone_baf_arm], ignore_index=True)
            df.append(temp_scdna_clone_baf)
    df = pd.concat(df, ignore_index=True)
    
    # save the output
    df.to_csv(argv.output, index=False)


if __name__ == '__main__':
    main()
