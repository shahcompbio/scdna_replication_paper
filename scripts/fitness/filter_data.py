import pandas as pd
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='unfiltered cn data for all cells')
    p.add_argument('output', type=str, help='cn data where all remaining cells have the same loci')

    return p.parse_args()


def filter_data(df):
    print("filtering based on gc, chrY")
    df2 = df.query('chr!="Y"')
    df2 = df2.query('gc > 0')

    mat2 = df2.pivot_table(index='cell_id', columns=['chr', 'start', 'end', 'gc'], values='reads')

    # drop loci that have >5% NaNs and then any remaining cells with >5% NaNs
    print("filtering bad loci")
    perc = 5.
    min_count = int(((100-perc)/100)*mat2.shape[0] + 1)
    mat3 = mat2.dropna(axis=1, thresh=min_count)

    print('filtering bad cells')
    perc = 5.
    min_count = int(((100-perc)/100)*mat3.shape[1] + 1)
    mat3 = mat3.dropna(axis=0, thresh=min_count)

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
    mat4 = df4.pivot_table(index='cell_id', columns=['chr', 'start', 'end', 'gc'], values='state')
    mat4 = mat4.fillna(axis=1, method='bfill').fillna(axis=1, method='ffill')
    df5 = mat4.reset_index().melt(id_vars='cell_id', value_name='state')
    print('merging state back into main df')
    df4.drop(columns=['state'], inplace=True)
    df6 = pd.merge(df4, df5, how='right')

    print('repeating process for copy')
    mat6 = df6.pivot_table(index='cell_id', columns=['chr', 'start', 'end', 'gc'], values='copy')
    mat6 = mat6.fillna(axis=1, method='bfill').fillna(axis=1, method='ffill')
    df7 = mat6.reset_index().melt(id_vars='cell_id', value_name='copy')
    print('merging copy back into main df')
    df6.drop(columns=['copy'], inplace=True)
    df8 = pd.merge(df6, df7, how='right')

    return df8


def merge_metrics(df, metrics):
    # subset to relevant metrics columns
    met = metrics[['cell_id', 'sample_id', 'library_id', 'quality', 'total_mapped_reads_hmmcopy', 'breakpoints', 'is_s_phase_prob']]

    # merge based on cell_id
    df = pd.merge(df, met, on='cell_id')

    # filter out any bad loci that get introduced from this merge
    # not sure why this happens but for some reason I need to do it
    mat = df.pivot_table(index='cell_id', columns=['chr', 'start', 'end', 'gc'], values='rpm')
    mat.dropna(axis=1, inplace=True)
    df2 = mat.reset_index().melt(id_vars='cell_id', value_name='rpm')
    df.drop(columns=['rpm'], inplace=True)
    df3 = pd.merge(df, df2, how='right')

    return df3



if __name__ == '__main__':
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')
    
    # make sure chromosome column is set to the appropriate dtype
    df['chr'] = df['chr'].astype(str)
    # for some reason it has to be set to a string when pivoting into a table
    # df['chr'] = df['chr'].astype('category')

    metrics = df[['cell_id', 'sample_id', 'library_id', 'quality', 'total_mapped_reads_hmmcopy', 'breakpoints', 'is_s_phase_prob']].drop_duplicates()

    # remove rpm column since I'll need to recalculate only using the good loci
    df = df.drop(columns=['rpm'])
    
    # remove appropriate cells and loci
    df = filter_data(df)

    # remove all columns that might contain NaNs
    df = df[['cell_id', 'chr', 'start', 'end', 'gc', 'reads', 'state', 'copy', 'rpm']]

    # merge long-form copy number info with cell metrics
    df = merge_metrics(df, metrics)

    df.to_csv(argv.output, sep='\t', index=False)
    