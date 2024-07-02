from argparse import ArgumentParser
import pandas as pd
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles
from scdna_replication_tools.assign_s_to_clones import assign_s_to_clones


def get_args():
    p = ArgumentParser()

    p.add_argument('cn', help='long-form cn_data for all cells in this dataset w/o clone assignments')
    p.add_argument('clones', help='clone assignments from the signatures paper')
    p.add_argument('dataset', help='name of this dataset')
    p.add_argument('assign_col', help='column to use for assigning S-phase cells to G1 clones')
    p.add_argument('cn_out', help='same as cn input but with clone_id brought in from tree or assigned via correlation for each cell')

    return p.parse_args()


def correct_normal_cells(df_s, df_g, thresh=0.90):
    ''' 
    Identify normal cells (>thresh% state==2 bins) and move them from df_s to df_g.
    This is necessary as normal cells are not included in the tree but they are not S-phase. 
    '''
    # compute the fraction of cn state==2 bins for each cell
    df_s['frac_state_2'] = df_s.groupby('cell_id')['state'].transform(lambda x: (x==2).sum() / len(x))
    # get cells with >thresh% state==2 bins and quality > 0.75
    # these are the normal cells which are confidently in G1/2-phase
    normal_cells = df_s.query('frac_state_2>@thresh & quality > 0.75')

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


def rename_none_clone(cn):
    """ 
    Rename clone with 'None' or 0 to the next available clone id. 
    For instance, if the clones are A, B, None, then rename None to C. 
    """
    # get all clone ids
    clones = cn['clone_id'].unique()
    # get the next available clone id not counting 'None'
    highest_clone = max([c for c in clones if c not in ['None', 0, '0']])
    # next clone id is the ascii character following the highest clone id
    next_clone = chr(ord(highest_clone) + 1)
    # rename 'None' clone to the next available clone id
    cn.loc[cn['clone_id'] == 'None', 'clone_id'] = next_clone
    cn.loc[cn['clone_id'] == 0, 'clone_id'] = next_clone
    cn.loc[cn['clone_id'] == '0', 'clone_id'] = next_clone
    return cn


def main():
    argv = get_args()

    cn = pd.read_csv(argv.cn)

    # load in clone assignments for each cell in this dataset
    # assignments are from most recent version of signatures paper analysis
    clones = pd.read_csv(argv.clones)
    clones = clones[['cell_id', 'clone_id']]

    # merge clone_id for the g1-phase cells (those with clone labels already)
    cn_g1 = pd.merge(cn, clones, on='cell_id')

    # rename 'None' clone to the next available clone id
    cn_g1 = rename_none_clone(cn_g1)
    
    # need to assign clone_ids for the S-phase cells who do not appear in cn_g1
    cn_s = cn.loc[~cn['cell_id'].isin(cn_g1['cell_id'].unique())]

    # correct normal cells and move them into a new clone within the G1/2-phase cells
    cn_s, cn_g1 = correct_normal_cells(cn_s, cn_g1)

    # compute conesensus clone profiles for assign_col
    clone_profiles = compute_consensus_clone_profiles(
        cn_g1, argv.assign_col, clone_col='clone_id', cell_col='cell_id', chr_col='chr',
        start_col='start', cn_state_col='state'
    )

    # assign S-phase cells to clones based on similarity of assign_col
    cn_s = assign_s_to_clones(
        cn_s, clone_profiles, col_name=argv.assign_col,
        clone_col='clone_id', cell_col='cell_id', chr_col='chr', start_col='start'
    )

    # note which cells were in the tree (had a clone assignment)
    cn_g1['in_tree'] = True
    cn_s['in_tree'] = False

    cn_out = pd.concat([cn_s, cn_g1], ignore_index=True)

    cn_out.to_csv(argv.cn_out, index=False)


if __name__ == '__main__':
    main()