import pandas as pd
import numpy as np
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells with relevent metrics columns')
    p.add_argument('permute_rate', type=float, help='fraction of G1/2-phase cells that should have their label swapped to S')
    p.add_argument('dataset', type=str, help='dataset id -- used to set a custom random seed for each dataset')
    p.add_argument('cn_output', type=str, help='cn data with the swapped labels')

    return p.parse_args()


def permute_cell_cycle(cn, rate):
    # copy over the true cell cycle states to new column before permuting
    cn['true_cell_cycle_state'] = cn['cell_cycle_state']

    # get a list of all the cells that are truly in G1/2-phase
    true_g_cells = cn.query('true_cell_cycle_state!="S"')['cell_id'].unique()

    # find number of cells to swap labels given the rate
    N = int(rate * len(true_g_cells))

    # get cell_ids for all the cells that need their labels swapped
    cells_to_swap = list(np.random.choice(true_g_cells, N, replace=False))

    # switch cell_cycle_state==S for all cells in the swap list
    cn.loc[cn['cell_id'].isin(cells_to_swap), 'cell_cycle_state'] = 'S'

    return cn



if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input, sep='\t')

    # set random seed corresponding to dataset id
    np.random.seed(ord(argv.dataset))

    # split by cell cycle
    cn = permute_cell_cycle(cn, argv.permute_rate)

    # return the permuted cn dataframe
    cn.to_csv(argv.cn_output, sep='\t', index=False)
    
