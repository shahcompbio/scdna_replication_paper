import pandas as pd
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='unfiltered cn data for all cells')
    p.add_argument('s_out', type=str, help='cn data where all remaining cells have the same loci -- S-phase')
    p.add_argument('g1_out', type=str, help='cn data where all remaining cells have the same loci -- G1-phase')
    p.add_argument('g2_out', type=str, help='cn data where all remaining cells have the same loci -- G2-phase')


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
    df4 = pd.merge(df2, df3, how='right')

    # convert to matrix and interpolate states
    mat4 = df4.pivot_table(index='cell_id', columns=['chr', 'start'], values='state')
    mat4 = mat4.fillna(axis=1, method='bfill').fillna(axis=1, method='ffill')
    df5 = mat4.reset_index().melt(id_vars='cell_id', value_name='state')
    df4.drop(columns=['state'], inplace=True)
    df6 = pd.merge(df4, df5, how='right')

    return df6


def subset_cell_cycle(df):
    df_s = df.query('cell_cycle_state=="S"')
    df_g1 = df.query('cell_cycle_state=="G1"')
    df_g2 = df.query('cell_cycle_state=="G2"')

    return df_s, df_g1, df_g2


def compute_rpm(df, input_col='reads', output_col='rpm', fill_columns=['cell_cycle_state', 'library_id']):
    df[output_col] = 0
    for cell_id, cell_cn in df.groupby('cell_id'):
        # compute reads per million for each cell
        cell_rpm = (cell_cn[input_col] * 1e6) / cell_cn[input_col].sum()
        df.loc[cell_cn.index, output_col] = cell_rpm

        # make sure all loci for this cell have same cell cycle state (no NaNs)
        # this is necessary because some bins only had NaNs filled for state & reads but not other columns
        for col in fill_columns:
            temp = cell_cn[col].fillna(method='bfill').fillna(method='ffill')
            df.loc[cell_cn.index, col] = temp

    return df


if __name__ == '__main__':
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')
    
    # remove appropriate cells and loci
    df = filter_data(df)

    # use raw read count to compute reads per million (each cell summs to 1 million reads)
    df = compute_rpm(df)

    # split based on flow sorted cell cycle state
    df_s, df_g1, df_g2 = subset_cell_cycle(df)

    df_s.to_csv(argv.s_out, sep='\t', index=False)
    df_g1.to_csv(argv.g1_out, sep='\t', index=False)
    df_g2.to_csv(argv.g2_out, sep='\t', index=False)
