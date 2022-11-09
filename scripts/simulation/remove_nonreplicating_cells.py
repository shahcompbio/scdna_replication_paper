import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('pyro_input', type=str, help='pyro model output for s-phase cells')
    p.add_argument('bulk_input', type=str, help='bulk model output for s-phase cells')
    p.add_argument('frac_rt_col', type=str, help='column name for the fraction of replicated bins in a cell')
    p.add_argument('pyro_rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('bulk_rep_col', type=str, help='column name for replicated status of each bin in bulk model')
    p.add_argument('pyro_output', type=str, help='pyro model for s-phase cells after filtering')
    p.add_argument('bulk_output', type=str, help='bulk model for s-phase cells after filtering')

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


def compute_ecf(cn_pyro, cn_bulk, frac_rt_col='cell_frac_rep'):
    ''' Find which cells have extreme cell fraction of replicated bins for both the pyro and bulk models. '''
    cell_metric_cols = [
        'cell_id', 'library_id', frac_rt_col, 
        'extreme_cell_frac', 'breakpoints', 
        'num_reads', 'madn', 'lrs', 
        'corrected_madn', 'corrected_breakpoints'
    ]
    cell_metrics_pyro = cn_pyro[cell_metric_cols].drop_duplicates()
    cell_metrics_bulk = cn_bulk[cell_metric_cols].drop_duplicates()

    ecf_pyro = cell_metrics_pyro[['cell_id', 'extreme_cell_frac']].drop_duplicates().reset_index(drop=True).rename(columns={'extreme_cell_frac': 'ecf_pyro'})
    ecf_bulk = cell_metrics_bulk[['cell_id', 'extreme_cell_frac']].drop_duplicates().reset_index(drop=True).rename(columns={'extreme_cell_frac': 'ecf_bulk'})
    ecf = pd.merge(ecf_pyro, ecf_bulk)

    return ecf


def remove_nonreplicating_cells(cn_pyro, cn_bulk, ecf, frac_rt_col='cell_frac_rep'):
    cell_metric_cols = [
        'cell_id', 'library_id', frac_rt_col, 
        'extreme_cell_frac', 'breakpoints', 
        'num_reads', 'madn', 'lrs', 
        'corrected_madn', 'corrected_breakpoints'
    ]
    cell_metrics_pyro = cn_pyro[cell_metric_cols].drop_duplicates()
    cell_metrics_bulk = cn_bulk[cell_metric_cols].drop_duplicates()

    # merge information about where ecf cells are located
    cn_pyro = pd.merge(cn_pyro, ecf)

    # use ecf status and cell cycle classifier features to nominate "bad" non-replicating cells
    # we are potentially removing cells with ecf in either version of running the pyro model
    extreme_cells = cn_pyro.loc[(cn_pyro['ecf_pyro']==True) | (cn_pyro['ecf_bulk']==True)]
    bad_cells_df = extreme_cells.loc[(extreme_cells['corrected_breakpoints']<0.0) | (extreme_cells['corrected_madn']<0.0)]
    bad_cells = bad_cells_df.cell_id.unique()

    cn_pyro_filtered = cn_pyro[~cn_pyro['cell_id'].isin(bad_cells)].reset_index(drop=True)
    cn_bulk_filtered = cn_bulk[~cn_bulk['cell_id'].isin(bad_cells)].reset_index(drop=True)

    return cn_pyro_filtered, cn_bulk_filtered


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn_pyro = pd.read_csv(argv.pyro_input, sep='\t')
    cn_bulk = pd.read_csv(argv.bulk_input, sep='\t')

    # compute the fraction of replicated bins within each cell
    cn_pyro = compute_cell_frac(cn_pyro, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.pyro_rep_col)
    cn_bulk = compute_cell_frac(cn_bulk, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.bulk_rep_col)

    # label which cells have extreme_cell_frac for both the pyro and bulk
    # versions of inferring scRT states
    ecf = compute_ecf(cn_pyro, cn_bulk, frac_rt_col=argv.frac_rt_col)

    # remove same set of "bad" cells from all 3 datasets if they seem to be nonreplicating
    cn_pyro_filtered, cn_bulk_filtered = remove_nonreplicating_cells(cn_pyro, cn_bulk, ecf, frac_rt_col=argv.frac_rt_col)

    # return one cn dataframe for each dataset
    cn_pyro_filtered.to_csv(argv.pyro_output, sep='\t', index=False)
    cn_bulk_filtered.to_csv(argv.bulk_output, sep='\t', index=False)
