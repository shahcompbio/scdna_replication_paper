from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from matplotlib.colors import ListedColormap
from scipy.optimize import curve_fit
from scgenome import cncluster
from matplotlib.patches import Patch
from scdna_replication_tools.calculate_twidth import compute_and_plot_twidth


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('pseduobulk', help='RT pseudobulk for this dataset')
    p.add_argument('dataset')
    p.add_argument('lamb', type=float, help='amount of sequencing noise')
    p.add_argument('A', type=float, help='steepness of inflection point when drawing RT state')
    p.add_argument('frac_rt_col', help='inferred fraction replicated for each cell')
    p.add_argument('true_frac_col', help='true fraction replicated for each cell')
    p.add_argument('rep_state', help='inferred replication state for each bin')
    p.add_argument('true_rep_state', help='true replication state for each bin')
    p.add_argument('infer_mode', help='pyro model or bulk')
    p.add_argument('output_tsv', help='table of all the computed t_width values')
    p.add_argument('output_curves', help='T-width curves of true and inferred rt_states')

    return p.parse_args()


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state'):
    ''' Compute the fraction of replicated bins for all cells in `cn` '''
    cn['extreme_cell_frac'] = False
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn[rep_state_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac
    return cn


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # compute fraction of replicated bins per cells 
    if argv.frac_rt_col not in df.columns:
        df = compute_cell_frac(df, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.rep_state)

    # merge in pseudobulk RT columns into df
    df2 = pd.read_csv(argv.pseduobulk, sep='\t')
    df = pd.merge(df, df2)

    # compute time from scheduled replication for each bin
    df['time_from_scheduled_rt'] = df['pseduobulk_hours'] - (df[argv.frac_rt_col] * 10.0)
    df['true_time_from_scheduled_rt'] = df['true_pseduobulk_hours'] - (df[argv.true_frac_col] * 10.0)

    # dataframe storing computed T-width values
    t_width_df = []

    fig, ax = plt.subplots(2, 2, figsize=(10, 10), tight_layout=True)
    ax = ax.flatten()

    title_second_line = 'dataset: {}, A: {}, lambda: {}, infer: {}'.format(
        argv.dataset, argv.A, argv.lamb, argv.infer_mode
    )

    # compute T-width with inferred scRT states
    ax[1], Tw = compute_and_plot_twidth(df, tfs_col='time_from_scheduled_rt', rs_col=argv.rep_state, title='Inferred scRT heterogeneity\n{}'.format(title_second_line), ax=ax[1])
    t_width_df.append(pd.DataFrame({
        'dataset': [argv.dataset], 'A': [argv.A], 'lambda': [argv.lamb], 'infer_mode': [argv.infer_mode], 'per_cell': [False], 'T-width': [Tw]
    }))

    # compute T-width for true scRT states
    ax[0], Tw = compute_and_plot_twidth(df, tfs_col='true_time_from_scheduled_rt', rs_col=argv.true_rep_state, title='True scRT heterogeneity\n{}'.format(title_second_line), ax=ax[0])
    t_width_df.append(pd.DataFrame({
        'dataset': [argv.dataset], 'A': [argv.A], 'lambda': [argv.lamb], 'infer_mode': ['ground_truth'], 'per_cell': [False], 'T-width': [Tw]
    }))

    # repeat the T-width calculations with per-cell==True
    # compute T-width with inferred scRT states
    ax[3], Tw = compute_and_plot_twidth(df, tfs_col='time_from_scheduled_rt', rs_col=argv.rep_state, title='Inferred scRT heterogeneity\n{}'.format(title_second_line), ax=ax[3], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({
        'dataset': [argv.dataset], 'A': [argv.A], 'lambda': [argv.lamb], 'infer_mode': [argv.infer_mode], 'per_cell': [True], 'T-width': [Tw]
    }))

    # compute T-width for true scRT states
    ax[2], Tw = compute_and_plot_twidth(df, tfs_col='true_time_from_scheduled_rt', rs_col=argv.true_rep_state, title='True scRT heterogeneity\n{}'.format(title_second_line), ax=ax[2], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({
        'dataset': [argv.dataset], 'A': [argv.A], 'lambda': [argv.lamb], 'infer_mode': ['ground_truth'], 'per_cell': [True], 'T-width': [Tw]
    }))

    # save figure of all the T-width curves
    fig.savefig(argv.output_curves, bbox_inches='tight')

    # save a table of all the computed T-width values
    t_width_df = pd.concat(t_width_df, ignore_index=True)
    t_width_df.to_csv(argv.output_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
