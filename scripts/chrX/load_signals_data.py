import pandas as pd
import numpy as np
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, help='chrX RT delays for sample profiles')
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


# load the scDNA signals results
def load_scdna(dataset):
    path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/schnapps-results/persample/{d}_hscn.csv.gz'.format(d=dataset)
    df = pd.read_csv(path, dtype={'chr': str})
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
        # always load the scDNA SIGNALS data
        temp_dna = load_scdna(d)
        # compute the BAF across the whole chromosome
        dna_baf_chr = compute_chrom_baf_df(temp_dna, modality='dna')
        # compute the BAF across each chromosome arm
        dna_baf_arm = compute_arm_baf_df(temp_dna, modality='dna')
        # concatenate the two dataframes
        temp_dna_baf = pd.concat([dna_baf_chr, dna_baf_arm], ignore_index=True)
        # try loading the scRNA SIGNALS data if available
        try:
            temp_rna = load_scrna(d)
            # compute the BAF across the whole chromosome
            rna_baf_chr = compute_chrom_baf_df(temp_rna, modality='rna')
            # compute the BAF across each chromosome arm
            rna_baf_arm = compute_arm_baf_df(temp_rna, modality='rna')
            # concatenate the two dataframes
            temp_rna_baf = pd.concat([rna_baf_chr, rna_baf_arm], ignore_index=True)
            temp_df = pd.merge(temp_rna_baf, temp_dna_baf)
        except:
            # fill scRNA columns with None if the dataset doesn't have scRNA data
            temp_df = temp_dna_baf.copy()
            temp_df['rna_baf_mean'] = None
            temp_df['rna_baf_std'] = None
            print('adding empty RNA columns for {}'.format(d))
        df.append(temp_df)
    df = pd.concat(df, ignore_index=True)
    

    # save the output
    df.to_csv(argv.output, index=False)


if __name__ == '__main__':
    main()
