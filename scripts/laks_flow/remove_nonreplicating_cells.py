import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('t47d_input', type=str, help='pyro model output for t47d s-phase cells')
    p.add_argument('GM18507_input', type=str, help='pyro model output for GM18507 s-phase cells')
    p.add_argument('joint_input', type=str, help='pyro model output for s-phase cells when run on both cell lines combined')
    p.add_argument('frac_rt_col', type=str, help='column denoting the fraction of replicated bins within each cell')
    p.add_argument('t47d_output', type=str, help='same as t47d input except with nonreplicating cells removed')
    p.add_argument('GM18507_output', type=str, help='same as GM18507 input except with nonreplicating cells removed')
    p.add_argument('joint_output', type=str, help='same as joint input except with nonreplicating cells removed')

    return p.parse_args()


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep'):
    ''' Compute the fraction of replicated bins for all cells in `cn`, noting which cells have extreme values. '''
    cn['extreme_cell_frac'] = False
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn['model_rep_state'].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac

        if temp_frac > 0.95 or temp_frac < 0.05:
            cn.loc[cell_cn.index, 'extreme_cell_frac'] = True
    return cn


def compute_ecf(cn_t, cn_gm, cn_all, frac_rt_col='cell_frac_rep'):
    ''' Find which cells have extreme cell fraction of replicated bins for both the joint model and the split by cell line model. '''
    cell_metric_cols = ['cell_id', 'sample_id', 'library_id', frac_rt_col, 'extreme_cell_frac',
                    'model_u', 'model_s_time', 'breakpoints', 'total_mapped_reads_hmmcopy',
                    'model_nb_r_s', 'model_a', 'rpm', 'madn', 'lrs',
                    'corrected_madn', 'corrected_breakpoints',]

    cell_metrics_all = cn_all[cell_metric_cols].drop_duplicates()
    cell_metrics_t = cn_t[cell_metric_cols].drop_duplicates()
    cell_metrics_gm = cn_gm[cell_metric_cols].drop_duplicates()

    ecf_all = cell_metrics_all[['cell_id', 'extreme_cell_frac']].drop_duplicates().reset_index(drop=True).rename(columns={'extreme_cell_frac': 'ecf_joint'})
    ecf_t = cell_metrics_t[['cell_id', 'extreme_cell_frac']].drop_duplicates().reset_index(drop=True).rename(columns={'extreme_cell_frac': 'ecf_split'})
    ecf_gm = cell_metrics_gm[['cell_id', 'extreme_cell_frac']].drop_duplicates().reset_index(drop=True).rename(columns={'extreme_cell_frac': 'ecf_split'})
    ecf = pd.merge(ecf_all, pd.concat([ecf_gm, ecf_t], ignore_index=True))

    return ecf


def remove_nonreplicating_cells(cn_t, cn_gm, cn_all, ecf, frac_rt_col='cell_frac_rep'):
    cn_split = pd.concat([cn_t, cn_gm], ignore_index=True)

    ccc_feature_cols = ['cell_id', 'madn', 'corrected_madn', 'lrs', 'cell_cycle_state', 
                    'library_id', 'extreme_cell_frac', frac_rt_col,
                    'total_mapped_reads_hmmcopy', 'breakpoints', 'corrected_breakpoints']
    cell_features_all = cn_all[ccc_feature_cols].drop_duplicates()
    cell_features_split = cn_split[ccc_feature_cols].drop_duplicates()

    # merge information about where ecf cells are located
    cn_all = pd.merge(cn_all, ecf)

    # use ecf status and cell cycle classifier features to nominate "bad" non-replicating cells
    # we are potentially removing cells with ecf in either version of running the pyro model
    extreme_cells = cn_all.loc[(cn_all['ecf_joint']==True) | (cn_all['ecf_split']==True)]
    bad_cells_df = extreme_cells.loc[(extreme_cells['corrected_breakpoints']<0.0) | (extreme_cells['corrected_madn']<0.0)]
    bad_cells = bad_cells_df.cell_id.unique()

    cn_all_filtered = cn_all[~cn_all['cell_id'].isin(bad_cells)].reset_index(drop=True)
    cn_t_filtered = cn_t[~cn_t['cell_id'].isin(bad_cells)].reset_index(drop=True)
    cn_gm_filtered = cn_gm[~cn_gm['cell_id'].isin(bad_cells)].reset_index(drop=True)

    return cn_all_filtered, cn_t_filtered, cn_gm_filtered


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn_t = pd.read_csv(argv.t47d_input, sep='\t')
    cn_gm = pd.read_csv(argv.GM18507_input, sep='\t')
    cn_all = pd.read_csv(argv.joint_input, sep='\t')

    # compute the fraction of replicated bins within each cell
    cn_all = compute_cell_frac(cn_all, frac_rt_col=argv.frac_rt_col)
    cn_t = compute_cell_frac(cn_t, frac_rt_col=argv.frac_rt_col)
    cn_gm = compute_cell_frac(cn_gm, frac_rt_col=argv.frac_rt_col)

    # label which cells have extreme_cell_frac for both the joint and split
    # versions of running the pyro model
    ecf = compute_ecf(cn_t, cn_gm, cn_all, frac_rt_col=argv.frac_rt_col)

    # remove same set of "bad" cells from all 3 datasets if they seem to be nonreplicating
    cn_all_filtered, cn_t_filtered, cn_gm_filtered = remove_nonreplicating_cells(cn_t, cn_gm, cn_all, ecf, frac_rt_col=argv.frac_rt_col)

    # return one cn dataframe for each dataset
    cn_t_filtered.to_csv(argv.t47d_output, sep='\t', index=False)
    cn_gm_filtered.to_csv(argv.GM18507_output, sep='\t', index=False)
    cn_all_filtered.to_csv(argv.joint_output, sep='\t', index=False)
