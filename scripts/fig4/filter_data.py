import pandas as pd
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='unfiltered cn data for all cells')
    p.add_argument('output', type=str, help='cn data where all remaining cells have the same loci')

    return p.parse_args()


def filter_data(df):
    print("filtering based on gc, ideal, chrY")
    df2 = df.query('ideal==True')
    df2 = df2.query('chr!="Y"')
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


def compute_rpm(df, input_col='reads', output_col='rpm'):
    df[output_col] = 0
    for cell_id, cell_cn in df.groupby('cell_id'):
        # compute reads per million for each cell
        cell_rpm = (cell_cn[input_col] * 1e6) / cell_cn[input_col].sum()
        df.loc[cell_cn.index, output_col] = cell_rpm

    return df


if __name__ == '__main__':
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')
    
    # remove appropriate cells and loci
    df = filter_data(df)

    # use raw read count to compute reads per million (each cell summs to 1 million reads)
    # df = compute_rpm(df)

    # remove all columns that might contain NaNs
    df = df[['cell_id', 'chr', 'start', 'end', 'gc', 'reads', 'state', 'copy', 'rpm']]

    df.to_csv(argv.output, sep='\t', index=False)
    