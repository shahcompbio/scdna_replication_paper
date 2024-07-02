import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', type=str, help='long-form csv results from PERT S-phase cells')
    p.add_argument('cn_g', type=str, help='long-form csv results from PERT G1/2-phase cells')
    p.add_argument('clone_col', type=str, help='column name of the clone ID in the input files')
    p.add_argument('output', help='table of cell cycle phase + clone counts')

    return p.parse_args()


def main():
    argv = get_args()

    # read in the data
    cn_s = pd.read_csv(argv.cn_s)
    cn_g = pd.read_csv(argv.cn_g)

    # columns to keep when dropping duplicates
    keep_cols = ['cell_id', argv.clone_col]

    s_cells = cn_s[keep_cols].drop_duplicates().reset_index(drop=True)
    g_cells = cn_g[keep_cols].drop_duplicates().reset_index(drop=True)

    # add cell cycle phase column
    s_cells['phase'] = 'S'
    g_cells['phase'] = 'G1/2'

    all_cells = pd.concat([s_cells, g_cells], ignore_index=True)

    # count the number of cells belonging to each clone_id + phase
    clone_counts = all_cells.groupby([argv.clone_col, 'phase']).size().reset_index(name='num_cells')


