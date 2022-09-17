from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles
from scdna_replication_tools.assign_s_to_clones import assign_s_to_clones


def get_args():
    p = ArgumentParser()

    p.add_argument('cn', help='long-form cn_data for all cells in this dataset w/o clone assignments')
    p.add_argument('clones', help='clone assignments from the fitness paper (applies to most datasets)')
    p.add_argument('clones_2295', help='clone assignments for just the 2295 dataset')
    p.add_argument('clones_SA1188', help='clone assignments for just the SA1188 dataset')
    p.add_argument('dataset', help='name of this dataset')
    p.add_argument('assign_col', help='column to use for assigning S-phase cells to G1 clones')
    p.add_argument('cn_out', help='same as cn input but with clone_id brought in from tree or assigned via correlation for each cell')

    return p.parse_args()


def main():
    argv = get_args()

    cn = pd.read_csv(argv.cn, sep='\t')

    if argv.dataset == '2295':
        clones = pd.read_csv(argv.clones_2295, sep='\t')
    elif argv.dataset == 'SA1188':
        clones = pd.read_csv(argv.clones_SA1188, sep='\t')
    else:
        # load in clones from fitness results
        clones = pd.read_csv(argv.clones)
        clones = clones.drop(columns=['V1', 'datatag', 'sample_id'])
        clones = clones.rename(columns={'single_cell_id': 'cell_id', 'letters': 'clone_id'})

        # remove the 'a' or 'b' suffix from SA906 cell IDs in the clone mapping file
        if 'SA906' in argv.dataset:
            clones['cell_id'] = clones['cell_id'].apply(lambda x: x.replace(argv.dataset, 'SA906'))

    # merge clone_id for the g1-phase cells (those with clone labels already)
    cn_g1 = pd.merge(cn, clones, on='cell_id')
    
    # need to assign clone_ids for the S-phase cells who do not appear in cn_g1
    # TODO: make sure this is right
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

    cn_out.to_csv(argv.cn_out, sep='\t', index=False)


if __name__ == '__main__':
    main()