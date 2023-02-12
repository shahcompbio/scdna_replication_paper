import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
import logging
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('cn_data', type=str, help='path to cnv_data.h5 file from cellrangercnv output')
    p.add_argument('clones', type=str, help='path to csv of hongyus treealign and clonealign clones')
    p.add_argument('output', type=str, help='long form tsv output for all loci and cells')
   
    return p.parse_args()


def counts_to_df(f, data_type='cnvs', col_name='state'):
    # loop through all chromosomes and concatenate the raw counts into a long-form pandas dataframe
    counts = []
    for chrom in f[data_type].keys():
        temp_counts = np.array(f[data_type][chrom][:])
        num_rows, num_bins = temp_counts.shape
        num_cells = int((num_rows + 1) / 2)  # there are 2N-1 rows for N cells
        # only take the first N rows for temp_counts
        temp_counts = temp_counts[:num_cells, :]
        
        # create a dataframe with the raw counts
        bin_id = [i for i in range(num_bins)] 
        cell_id = [str(i) for i in range(num_cells)]
        temp_df = pd.DataFrame(temp_counts, columns=bin_id, index=cell_id)
        
        # melt the dataframe
        temp_df = temp_df.reset_index().melt(id_vars='index', var_name='bin_id', value_name=col_name)
        temp_df['chr'] = chrom.replace('chr', '')

        # rename the index column to cell_id
        temp_df.rename(columns={'index': 'cell_id'}, inplace=True)

        counts.append(temp_df)
    counts = pd.concat(counts, ignore_index=True)

    return counts


def add_genome_tracks(f):
    df = []
    # loop through all chromosomes and concatenate the raw counts into a long-form pandas dataframe
    for chrom in f['cnvs'].keys():
        # find the number of bins in this chromosome using the cnvs matrix
        temp_counts = np.array(f['cnvs'][chrom][:])
        _, num_bins = temp_counts.shape
        bin_id = [i for i in range(num_bins)]

        # create a temporaray dataframe for this chromosome
        temp_df = pd.DataFrame()
        temp_df['bin_id'] = bin_id
        temp_df['chr'] = chrom.replace('chr', '')
        temp_df['gc'] = np.array(f['genome_tracks']['gc_fraction'][chrom])
        temp_df['is_mappable'] = np.array(f['genome_tracks']['is_mappable'][chrom])
        temp_df['n_fraction'] = np.array(f['genome_tracks']['n_fraction'][chrom])

        # add the start and end positions of each bin using the bin size
        bin_size = int(np.array(f['constants']['bin_size']))
        temp_df['start'] = temp_df['bin_id'].astype(int) * bin_size + 1
        temp_df['end'] = temp_df['start'] + bin_size - 1
    
        # append the temporary dataframe to the main dataframe
        df.append(temp_df)
    # concatenate all the chromsome dataframes into one dataframe
    df = pd.concat(df, ignore_index=True)
    return df


def parse_h5_input(f):
    '''
    This function takes in the h5 file from cellrangercnv and returns a long-form dataframe with the following columns:
    cell_id, bin_id, chr, start, end, gc, state, reads, reads_norm
    '''
    print('Parsing h5 file')
    # load the cnvs
    cnvs = counts_to_df(f, data_type='cnvs', col_name='state')
    print('Loaded cnvs')
    # load gc, start, end, is_mappable, n_fraction for each locus
    genome_tracks = add_genome_tracks(f)
    print('Loaded genome tracks')
    # merge the cnvs and genome_tracks dataframes
    cnvs = cnvs.merge(genome_tracks)
    print('Merged cnvs and genome tracks')
    # load raw_counts
    raw_counts = counts_to_df(f, data_type='raw_counts', col_name='reads')
    print('Loaded raw counts')
    # load normalized counts
    norm_counts = counts_to_df(f, data_type='normalized_counts', col_name='reads_norm')
    print('Loaded normalized counts')
    # merge the raw counts, normalized counts and cnvs dataframes
    df = pd.merge(raw_counts, norm_counts)
    df = pd.merge(df, cnvs)
    print('Merged raw counts, normalized counts and cnvs')
    # filter out the unmappable regions
    df = df[(df['is_mappable'] == 1)]
    print('Filtered out unmappable regions')
    # convert df['cell_id'] to int64
    df['cell_id'] = df['cell_id'].astype(int)
    # drop n_fraction and is_mappable columns as we've already filtered out the unmappable bins
    df.drop(columns=['n_fraction', 'is_mappable'], inplace=True)
    print('Dropped n_fraction and is_mappable columns')

    return df


def parse_cell_metrics(f):
    '''
    This function takes in the h5 file from cellrangercnv and returns a dataframe precomputed per-cell metrics
    '''
    # loop through all the per cell summary metrics and add them as columns to a pandas dataframe
    metrics_df = pd.DataFrame()
    for metric in f['per_cell_summary_metrics'].keys():
        temp_metric = np.array(f['per_cell_summary_metrics'][metric])
        metrics_df[metric] = temp_metric

    # convert the barcode column to a string without the b'*' format
    metrics_df['barcode'] = metrics_df['barcode'].str.decode('utf-8')

    return metrics_df


def add_clones_to_metrics(argv, metrics_df):
    '''
    This function takes in the metrics dataframe and adds the clone information to it
    '''
    # load the clones dataframe
    clones = pd.read_csv(argv.clones)
    # drop the first column which is just the index in the csv file
    clones.drop(columns=['Unnamed: 0'], inplace=True)
    # rename the cell_id column to barcode
    clones.rename(columns={'cell_id': 'barcode'}, inplace=True)
    # rename clonealign_clone_id to clone_id
    clones.rename(columns={'clonealign_clone_id': 'clone_id'}, inplace=True)
    # if clone_id column is of dtype int, add a 'clone_id' column that converts 1->A, 2->B, 3->C, etc.
    if clones['clone_id'].dtype == 'int64':
        clones['clone_id'] = clones['clone_id'].apply(lambda x: chr(x + 64))

    # merge the clones dataframe with the metrics dataframe
    # keep rows with no clone assignment (NaN)
    metrics_df = pd.merge(metrics_df, clones, how='left')
    # fill in the missing clone assignments with -1 to show it is outside of the tree
    metrics_df['clonealign_tree_id'].fillna('node_-1', inplace=True)
    metrics_df['clone_id'].fillna('not_in_tree', inplace=True)

    return metrics_df


def remove_state0_cells(df):
    '''
    This function takes in the long-form dataframe and removes cells with more than 7500 of their bins in state 0
    '''
    # count the number of state==0 bins in each cell
    # use this feature to remove cells with too many state==0 bins from the set of S and low quality cells
    df['num_state_0_bins'] = df.groupby('cell_id')['state'].transform(lambda x: (x == 0).sum())

    # remove cells with more than 7500 state==0 bins
    df_out = df[df['num_state_0_bins'] < 7500]

    return df_out


def filter_high_cn_loci(df, thresh=12):
    '''
    This function takes in the long-form dataframe and filters out the loci that have state>12 bins in the population of high quality G1/2-phase cells
    '''
    # pivot to a wide-form dataframe of cells x bins
    mat = df.pivot_table(index='cell_id', columns=['chr', 'start'], values='state')

    # remove any columns from the matrix that have values >thresh
    mat2 = mat.loc[:, mat.max() <= thresh]

    # melt back to long form and merge with the original dataframe to recover the other columns
    df2 = mat2.reset_index().melt(id_vars='cell_id', var_name=['chr', 'start'], value_name='state')
    df_out = df2.merge(df, how='left')

    return df_out


def filter_other_bad_cells_and_loci(df):
    print("filtering based on gc, chrY")
    df2 = df.query('chr!="Y"')
    df2 = df2.query('gc > 0')

    mat2 = df2.pivot_table(index='cell_id', columns=['chr', 'start'], values='reads')

    # drop cells that have >5% NaNs and then any remaining loci with >5% NaNs
    # NaNs in this matrix are caused by loci being present in some cells but not others
    print("filtering bad cells")
    perc = 5.
    min_count = int(((100-perc)/100)*mat2.shape[1] + 1)
    mat3 = mat2.dropna(axis=0, thresh=min_count)

    print('filtering bad loci')
    perc = 5.
    min_count = int(((100-perc)/100)*mat3.shape[0] + 1)
    mat3 = mat3.dropna(axis=1, thresh=min_count)

    # fill the remaining NaNs with the neighboring loci for that same cell
    print('filling remaining NaNs')
    mat3 = mat3.fillna(axis=1, method='bfill').fillna(axis=1, method='ffill')
    # compute reads per million for each cell
    print('computing rpm')
    mat_rpm = mat3.div(mat3.sum(axis=1), axis=0) * 1e6

    # limit original df to just include melted loci
    print('melting and merging reads and rpm')
    df3 = mat3.reset_index().melt(id_vars='cell_id', value_name='reads')
    df_rpm = mat_rpm.reset_index().melt(id_vars='cell_id', value_name='rpm')

    # combine filtered reads and rpm into one df
    df3 = pd.merge(df3, df_rpm)

    # merge using the interpolated read counts from df3
    print('merging reads and rpm back into main df')
    df2.drop(columns=['reads'], inplace=True)
    df4 = pd.merge(df2, df3, how='right')

    # convert to matrix and interpolate states
    print('repeating process for state')
    mat4 = df4.pivot_table(index='cell_id', columns=['chr', 'start'], values='state')
    mat4 = mat4.fillna(axis=1, method='bfill').fillna(axis=1, method='ffill')
    df5 = mat4.reset_index().melt(id_vars='cell_id', value_name='state')
    print('merging state back into main df')
    df4.drop(columns=['state'], inplace=True)
    df6 = pd.merge(df4, df5, how='right')

    # merge with the original dataframe to recover the cell metric columns
    df_out = df6.merge(df, how='left')

    return df_out


def main():
    argv = get_args()
    print('Parsed the arguments')

    # read in the cnv data from the h5 file
    f = h5py.File(argv.cn_data, 'r')
    print('loaded the h5 file')

    # parse the h5 file into a long-form cn dataframe
    df = parse_h5_input(f)
    print('Parsed the h5 file into a long-form dataframe')

    # parse the h5 file into a dataframe with per-cell metrics
    metrics_df = parse_cell_metrics(f)
    print('Parsed the h5 file into a dataframe with per-cell metrics')

    # close the h5 file
    f.close()

    # add the clone information to the metrics dataframe
    metrics_df = add_clones_to_metrics(argv, metrics_df)
    print('Added the clone information to the metrics dataframe')

    # merge metrics_df with df
    df2 = pd.merge(df, metrics_df)
    print('Merged the metrics dataframe with the long-form dataframe')

    # remove cells with more than 7500 state==0 bins
    df3 = remove_state0_cells(df2)
    print('Removed cells with more than 7500 state==0 bins')

    # filter out the loci that have state>12 bins in the population of G1/2-phase cells
    df4 = filter_high_cn_loci(df3)
    print('Filtered out the loci that have state>12 bins in the population of G1/2-phase cells')

    # filter other bad cells and loci
    df5 = filter_other_bad_cells_and_loci(df4)

    # save the dataframe to a csv file
    df5.to_csv(argv.output, index=False)


if __name__ == '__main__':
    main()
