import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    
    parser.add_argument('cn_s', type=str, help='path to long-form PERT output of S-phase cells')
    parser.add_argument('cn_g', type=str, help='path to long-form PERT output of G1/2-phase cells')
    parser.add_argument('clone_counts', type=str, help='csv file with cell cycle clone counts')
    parser.add_argument('sample_counts', type=str, help='csv file with cell cycle counts across the whole sample')
    
    return parser.parse_args()


def main():
    argv = get_args()

    # read in the long-form PERT output for S-phase and G1/2-phase cells
    cn_s = pd.read_csv(argv.cn_s)
    cn_g = pd.read_csv(argv.cn_g)

    # reduce to just the essential per-cell columns
    clone_col = 'assigned_clone_id'
    keep_cols = ['cell_id', clone_col]
    s_cells = cn_s[keep_cols].drop_duplicates().reset_index(drop=True)
    g_cells = cn_g[keep_cols].drop_duplicates().reset_index(drop=True)

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
    s_counts = s_cells.groupby([clone_col]).size().reset_index(name='dna_num_cells_s')
    g_counts = g_cells.groupby([clone_col]).size().reset_index(name='dna_num_cells_g1')

    # merge the two dataframes giving the phase counts for each clone
    clone_counts = s_counts.merge(g_counts, on=clone_col, how='outer')

    # for each clone (row), calculate the fraction of cells in S-phase and G1 phase
    clone_counts['dna_frac_s'] = clone_counts['dna_num_cells_s'] / (clone_counts['dna_num_cells_s'] + clone_counts['dna_num_cells_g1'])
    clone_counts['dna_frac_g1'] = clone_counts['dna_num_cells_g1'] / (clone_counts['dna_num_cells_s'] + clone_counts['dna_num_cells_g1'])

    # save the output clone and sample counts as csv files
    clone_counts.to_csv(argv.clone_counts, index=False)
    sample_counts.to_csv(argv.sample_counts, index=False)



if __name__ == '__main__':
    main()
