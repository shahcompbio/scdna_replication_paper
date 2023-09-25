import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells with relevent metrics columns')
    p.add_argument('s_out', type=str, help='cn data for S-phase cells')
    p.add_argument('g_out', type=str, help='cn data for G1/2-phase cells')

    return p.parse_args()


def remove_normal_cells(df_s, df_g):
    ''' 
    Identify normal cells (>95% state==2 bins) and move them from df_s to df_g.
    This is necessary as normal cells are not included in the tree but they are not S-phase. 
    '''
    # compute the fraction of cn state==2 bins for each cell
    df_s['frac_state_2'] = df_s.groupby('cell_id')['state'].transform(lambda x: (x==2).sum() / len(x))
    # get cells with >95% state==2 bins
    normal_cells = df_s.query('frac_state_2>0.95')

    # remove normal cells from df_s
    df_s = df_s.loc[~df_s['cell_id'].isin(normal_cells['cell_id'])].reset_index(drop=True)

    # find the next available clone id for normal cells
    # get all clone ids
    clones = df_g['clone_id'].unique()
    # get the next available clone id not counting 'None'
    highest_clone = max([c for c in clones if c != 'None'])
    # next clone id is the ascii character following the highest clone id
    next_clone = chr(ord(highest_clone) + 1)
    # rename the normal cells to the next available clone id
    normal_cells['clone_id'] = next_clone

    # concatenate normal_cells to df_g
    df_g = pd.concat([df_g, normal_cells], axis=0).reset_index(drop=True)

    return df_s, df_g


def split_cell_cycle(cn):
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
    cn = pd.read_csv(argv.cn_input)

    # split by cell cycle
    df_s, df_g = split_cell_cycle(cn)

    # remove normal cells
    df_s, df_g = remove_normal_cells(df_s, df_g)

    # return one cn dataframe for each cell cycle state
    df_s.to_csv(argv.s_out, index=False)
    df_g.to_csv(argv.g_out, index=False)

