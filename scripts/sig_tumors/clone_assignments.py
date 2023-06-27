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


def rename_none_clone(cn):
    """ 
    Rename clone with 'None' to the next available clone id. 
    For instance, if the clones are A, B, None, then rename None to C. 
    """
    # get all clone ids
    clones = cn['clone_id'].unique()
    # get the next available clone id not counting 'None'
    highest_clone = max([c for c in clones if c != 'None'])
    # next clone id is the ascii character following the highest clone id
    next_clone = chr(ord(highest_clone) + 1)
    # rename 'None' clone to the next available clone id
    cn.loc[cn['clone_id'] == 'None', 'clone_id'] = next_clone
    return cn


def main():
    argv = get_args()

    cn = pd.read_csv(argv.cn)

    # load in clone assignments for each cell in this dataset
    # assignments are from most recent version of signatures paper analysis
    clones = pd.read_csv(argv.clones, sep='\t')

    # merge clone_id for the g1-phase cells (those with clone labels already)
    cn_g1 = pd.merge(cn, clones, on='cell_id')

    # rename 'None' clone to the next available clone id
    cn_g1 = rename_none_clone(cn_g1)
    
    # need to assign clone_ids for the S-phase cells who do not appear in cn_g1
    cn_s = cn.loc[~cn['cell_id'].isin(cn_g1['cell_id'].unique())]

    # note which cells were in the tree (had a clone assignment)
    cn_g1['in_tree'] = True
    cn_s['in_tree'] = False

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

    cn_out = pd.concat([cn_s, cn_g1], ignore_index=True)

    cn_out.to_csv(argv.cn_out, index=False)


if __name__ == '__main__':
    main()