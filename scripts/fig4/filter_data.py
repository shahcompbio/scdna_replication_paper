import pandas as pd
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='unfiltered cn data for all cells')
    p.add_argument('output', type=str, help='cn data where all remaining cells have the same loci')

    return p.parse_args()


def filter_data(df):
    df2 = df.query('ideal==True')
    df2 = df2.query('chr!="Y"')
    df2 = df2.query('gc > 0')

    mat2 = df2.pivot_table(index='cell_id', columns=['chr', 'start'], values='reads')

    # drop loci that have >5% NaNs and then any remaining cells with >5% NaNs
    perc = 5.
    min_count = int(((100-perc)/100)*mat2.shape[0] + 1)
    mat3 = mat2.dropna(axis=1, thresh=min_count)

    perc = 5.
    min_count = int(((100-perc)/100)*mat3.shape[1] + 1)
    mat3 = mat3.dropna(axis=0, thresh=min_count)

    # fill the remaining NaNs with the neighboring loci for that same cell
    mat3 = mat3.fillna(axis=1, method='bfill').fillna(axis=1, method='ffill')

    # limit original df to just include melted loci
    df3 = mat3.reset_index().melt(id_vars='cell_id', value_name='reads')

    # merge using the interpolated read counts from df3
    df2.drop(columns=['reads'], inplace=True)
    df4 = pd.merge(df2, df3)

    # convert to matrix and interpolate states
    mat4 = df4.pivot_table(index='cell_id', columns=['chr', 'start'], values='state')
    mat4 = mat4.fillna(axis=1, method='bfill').fillna(axis=1, method='ffill')
    df5 = mat4.reset_index().melt(id_vars='cell_id', value_name='state')
    df4.drop(columns=['state'], inplace=True)
    df6 = pd.merge(df4, df5)

    return df6



if __name__ == '__main__':
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')
    df = filter_data(df)
    df.to_csv(argv.output, sep='\t', index=False)
