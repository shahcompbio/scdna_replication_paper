import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='pyro model output for s-phase cells')
    p.add_argument('frac_rt_col', type=str, help='column name for the fraction of replicated bins in a cell')
    p.add_argument('rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('output', type=str, help='pyro model for s-phase cells after filtering')

    return p.parse_args()


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state'):
    ''' Compute the fraction of replicated bins for all cells in `cn`, noting which cells have extreme values. '''
    cn['extreme_cell_frac'] = False
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn[rep_state_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac

        if temp_frac > 0.95 or temp_frac < 0.05:
            cn.loc[cell_cn.index, 'extreme_cell_frac'] = True
    return cn


def remove_nonreplicating_cells(cn, frac_rt_col='cell_frac_rep'):
    # use extreme_cell_frac status and cell cycle classifier features to nominate "bad" non-replicating cells
    extreme_cells = cn.loc[(cn['extreme_cell_frac']==True)]
    bad_cells_df = extreme_cells.loc[(extreme_cells['corrected_breakpoints']<0.0) | (extreme_cells['corrected_madn']<0.0)]
    bad_cells = bad_cells_df.cell_id.unique()

    cn_filtered = cn[~cn['cell_id'].isin(bad_cells)].reset_index(drop=True)

    return cn_filtered


if __name__ == '__main__':
    argv = get_args()

    # load in unfiltered S-phase cells
    cn = pd.read_csv(argv.input, sep='\t')

    # compute the fraction of replicated bins within each cell
    cn = compute_cell_frac(cn, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.rep_col)

    # remove same set of "bad" cells from all 3 datasets if they seem to be nonreplicating
    cn_filtered = remove_nonreplicating_cells(cn, frac_rt_col=argv.frac_rt_col)

    # return the filtered set of S-phase cells
    cn_filtered.to_csv(argv.output, sep='\t', index=False)
