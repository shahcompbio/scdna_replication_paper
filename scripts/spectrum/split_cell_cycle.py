import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells with relevent metrics columns')
    p.add_argument('s_out', type=str, help='cn data for S-phase cells')
    p.add_argument('g_out', type=str, help='cn data for G1/2-phase cells')

    return p.parse_args()


def correct_normal_cells(df_s, df_g, thresh=0.90):
    ''' 
    Identify normal cells (>thresh% state==2 bins) and move them from df_s to df_g.
    This is necessary as normal cells are not included in the tree but they are not S-phase. 
    '''
    # compute the fraction of cn state==2 bins for each cell
    df_s['frac_state_2'] = df_s.groupby('cell_id')['state'].transform(lambda x: (x==2).sum() / len(x))
    # get cells with >thresh% state==2 bins, quality > 0.75, corrected_madn < 0.05, and corrected_breakpoints < 0
    # these are the normal cells which are confidently in G1/2-phase
    normal_cells = df_s.query('frac_state_2>@thresh & quality > 0.75 & corrected_madn < 0.05 & corrected_breakpoints < 0')

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

    # change all the other normal cells initialized as S-phase to have the same clone_id
    # these are the cells that have frac_state_2 but are still in df_s
    df_s.loc[df_s['frac_state_2'] > thresh, 'clone_id'] = next_clone

    # concatenate normal_cells to df_g
    df_g = pd.concat([df_g, normal_cells], axis=0).reset_index(drop=True)

    return df_s, df_g


def correct_high_gc_errors(df_s, df_g):
    ''' 
    Some cells might be excluded from the tree due to high GC content which would cause the Laks classifier
    to misclassify them as S-phase but shouldn't show up in the other metrics. This function will identify those
    cells and move them from df_s to df_g. 
    '''
    # subset to the cells in df_s with quality > 0.75, and is_s_phase == True
    high_gc_cells = df_s.query('quality>0.75 & is_s_phase==True')

    # remove the high_gc_cells from df_s
    df_s = df_s.loc[~df_s['cell_id'].isin(high_gc_cells['cell_id'])].reset_index(drop=True)

    # add the high_gc_cells to df_g
    df_g = pd.concat([df_g, high_gc_cells], axis=0).reset_index(drop=True)

    return df_s, df_g


def split_by_tree(cn):
    # make sure column dtypes are correct
    cn['chr'] = cn['chr'].astype(str)
    cn['start'] = cn['start'].astype(int)
    cn['end'] = cn['end'].astype(int)
    cn['cell_id'] = cn['cell_id'].astype(str)

    # split into initial cell cycle phases based on whether the cell is in the tree or not
    # this is good start because the tree consists of high-confidence G1/2-phase tumor cells
    df_s = cn.query('in_tree==False').reset_index(drop=True)
    df_g = cn.query('in_tree==True').reset_index(drop=True)

    return df_s, df_g


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input)

    # start by whether the cell is in the tree or not
    df_s, df_g = split_by_tree(cn)

    # # correct for normal cells
    # df_s, df_g = correct_normal_cells(df_s, df_g)

    # correct for high GC bias errors
    df_s, df_g = correct_high_gc_errors(df_s, df_g)

    # return one cn dataframe for each cell cycle state
    df_s.to_csv(argv.s_out, index=False)
    df_g.to_csv(argv.g_out, index=False)

