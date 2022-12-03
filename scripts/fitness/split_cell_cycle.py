import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells with relevent metrics columns')
    p.add_argument('dataset')
    p.add_argument('s_out', type=str, help='cn data for S-phase cells')
    p.add_argument('g_out', type=str, help='cn data for G1/2-phase cells')

    return p.parse_args()


def split_cell_cycle(cn):
    ## TODO: come up with scheme more complex than just in/out tree
    ## maybe the cells in the top quartile of ccc features within each library?
    # make sure column dtypes are correct
    cn['chr'] = cn['chr'].astype(str)
    cn['start'] = cn['start'].astype(int)
    cn['end'] = cn['end'].astype(int)
    cn['cell_id'] = cn['cell_id'].astype(str)

    cn_g = cn.loc[(cn['quality']>0.9) & (cn['corrected_breakpoints']<0)].reset_index(drop=True)
    cn_s = cn.loc[~cn['cell_id'].isin(cn_g['cell_id'].unique())].reset_index(drop=True)

    return cn_s, cn_g


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input, sep='\t')

    # split by cell cycle
    cn_s, cn_g = split_cell_cycle(cn)

    # return one cn dataframe for each cell cycle state
    cn_s.to_csv(argv.s_out, sep='\t', index=False)
    cn_g.to_csv(argv.g_out, sep='\t', index=False)

