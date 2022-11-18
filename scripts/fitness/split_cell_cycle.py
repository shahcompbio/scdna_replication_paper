import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells with relevent metrics columns')
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

    # subset cn based on cell cycle state
    df_s = cn.query('in_tree==False').reset_index(drop=True)
    df_g = cn.query('in_tree==True').reset_index(drop=True)

    return df_s, df_g



if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input, sep='\t')

    # split by cell cycle
    df_s, df_g = split_cell_cycle(cn)

    ## TODO: recluster the G1/2-phase cells (new clone_id's)


    ## TODO: assign the S-phase cells to these new clusters

    

    # return one cn dataframe for each cell cycle state
    df_s.to_csv(argv.s_out, sep='\t', index=False)
    df_g.to_csv(argv.g_out, sep='\t', index=False)

