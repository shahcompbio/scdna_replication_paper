import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    
    parser.add_argument('cn_s', type=str, help='path to long-form PERT output of S-phase cells')
    parser.add_argument('cn_g', type=str, help='path to long-form PERT output of G1/2-phase cells')
    parser.add_argument('rna_cells', type=str, help='csv file with cell_id, phase, and clone_id for scRNA cells')
    parser.add_argument('clone_counts', type=str, help='csv file with cell cycle clone counts')
    parser.add_argument('sample_counts', type=str, help='csv file with cell cycle counts across the whole sample')
    
    return parser.parse_args()


def compute_rna_clone_counts(rna_cells, clone_col='clone_id'):
    """ 
    Provided a dataframe with cell_id, phase, and clone_id for scRNA cells,
    return a dataframe with the number of cells in each phase for each clone_id
    """
    # rename the phase and clone_id columns to to be consistent with DNA data
    rna_cells.rename(columns={'Phase': 'phase', 'clone_id': clone_col}, inplace=True)

    # split the rna_cells dataframe into three dataframes, one for each phase
    rna_cells_s = rna_cells[rna_cells['phase'] == 'S']
    rna_cells_g1 = rna_cells[rna_cells['phase'] == 'G1']
    rna_cells_g2m = rna_cells[rna_cells['phase'] == 'G2M']

    # count the number of rna cells in each phase
    s_rna_counts = rna_cells_s.groupby([clone_col]).size().reset_index(name='rna_num_cells_s')
    g1_rna_counts = rna_cells_g1.groupby([clone_col]).size().reset_index(name='rna_num_cells_g1')
    g2m_rna_counts = rna_cells_g2m.groupby([clone_col]).size().reset_index(name='rna_num_cells_g2m')

    # merge the three dataframes giving the phase counts for each clone
    rna_clone_counts = s_rna_counts.merge(g1_rna_counts, on=clone_col, how='outer')
    rna_clone_counts = rna_clone_counts.merge(g2m_rna_counts, on=clone_col, how='outer')

    # for each clone (row) in the RNA counts, calculate the fraction of cells belonging to each phase
    rna_clone_counts['rna_frac_s'] = rna_clone_counts['rna_num_cells_s'] / (rna_clone_counts['rna_num_cells_s'] + rna_clone_counts['rna_num_cells_g1'] + rna_clone_counts['rna_num_cells_g2m'])
    rna_clone_counts['rna_frac_g1'] = rna_clone_counts['rna_num_cells_g1'] / (rna_clone_counts['rna_num_cells_s'] + rna_clone_counts['rna_num_cells_g1'] + rna_clone_counts['rna_num_cells_g2m'])
    rna_clone_counts['rna_frac_g2m'] = rna_clone_counts['rna_num_cells_g2m'] / (rna_clone_counts['rna_num_cells_s'] + rna_clone_counts['rna_num_cells_g1'] + rna_clone_counts['rna_num_cells_g2m'])

    return rna_clone_counts


def main():
    argv = get_args()

    # read in the long-form PERT output for S-phase and G1/2-phase cells and scRNA cell phase data
    cn_s = pd.read_csv(argv.cn_s)
    cn_g = pd.read_csv(argv.cn_g)
    rna_cells = pd.read_csv(argv.rna_cells)

    clone_col = 'assigned_clone_id'
    treealign_clone_col = 'clonealign_tree_id'

    # find the proper mapping from 'clone_id' to 'clonealign_tree_id' 
    tree_to_clone_map = cn_g[['clone_id', treealign_clone_col]].drop_duplicates().reset_index(drop=True)
    # rename the 'clone_id' column to clone_col
    tree_to_clone_map.rename(columns={'clone_id': clone_col}, inplace=True)

    # reduce to just the essential per-cell columns
    keep_cols = ['cell_id', clone_col]
    s_cells = cn_s[keep_cols].drop_duplicates().reset_index(drop=True)
    g_cells = cn_g[keep_cols].drop_duplicates().reset_index(drop=True)

    print('s_cells:', s_cells.head(), sep='\n')
    print('g_cells:', g_cells.head(), sep='\n')
    print('tree_to_clone_map:', tree_to_clone_map.head(), sep='\n')

    # merge the tree align clone ID column into each dataframe
    s_cells = s_cells.merge(tree_to_clone_map, on=clone_col, how='left')
    g_cells = g_cells.merge(tree_to_clone_map, on=clone_col, how='left')

    print('s_cells:', s_cells.head(), sep='\n')
    print('g_cells:', g_cells.head(), sep='\n')

    # count the number of cells in each phase for the whole sample
    sample_counts = pd.DataFrame(
        {
            'dna_num_cells_s': [len(s_cells)],
            'dna_num_cells_g1': [len(g_cells)],
            'dna_frac_s': [len(s_cells) / (len(s_cells) + len(g_cells))],
            'dna_frac_g1': [len(g_cells) / (len(s_cells) + len(g_cells))],
        }
    )

    # count the number of cells belonging to each clone_id within each phase
    s_counts = s_cells.groupby([treealign_clone_col]).size().reset_index(name='dna_num_cells_s')
    g_counts = g_cells.groupby([treealign_clone_col]).size().reset_index(name='dna_num_cells_g1')

    print('s_counts:', s_counts.head(), sep='\n')
    print('g_counts:', g_counts.head(), sep='\n')

    # merge the two dataframes giving the phase counts for each clone
    clone_counts = s_counts.merge(g_counts, on=treealign_clone_col, how='outer')

    # for each clone (row), calculate the fraction of cells in S-phase and G1 phase
    clone_counts['dna_frac_s'] = clone_counts['dna_num_cells_s'] / (clone_counts['dna_num_cells_s'] + clone_counts['dna_num_cells_g1'])
    clone_counts['dna_frac_g1'] = clone_counts['dna_num_cells_g1'] / (clone_counts['dna_num_cells_s'] + clone_counts['dna_num_cells_g1'])

    print('clone_counts:', clone_counts.head(), sep='\n')

    # compute the number of cells in each phase for each clone_id in the scRNA data
    rna_clone_counts = compute_rna_clone_counts(rna_cells, clone_col=treealign_clone_col)

    print('rna_clone_counts:', rna_clone_counts.head(), sep='\n')

    # merge the RNA and DNA clone counts
    clone_counts = clone_counts.merge(rna_clone_counts, on=treealign_clone_col, how='inner')

    print('clone_counts:', clone_counts.head(), sep='\n')

    # save the output clone and sample counts as csv files
    clone_counts.to_csv(argv.clone_counts, index=False)
    sample_counts.to_csv(argv.sample_counts, index=False)



if __name__ == '__main__':
    main()
