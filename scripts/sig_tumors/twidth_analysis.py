from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scdna_replication_tools.calculate_twidth import compute_and_plot_twidth


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('bulk_rt', help='pseudobulk RT information')
    p.add_argument('dataset')
    p.add_argument('frac_rt_col', help='column denoting the fraction of replicated loci per cell (its time in S-phase)')
    p.add_argument('rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('output_csv', help='table of all the computed t_width values')
    p.add_argument('output_pdf', help='T-width curve using inferred scRT')

    return p.parse_args()


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s)

    bulk_rt = pd.read_csv(argv.bulk_rt)
    bulk_rt = bulk_rt[['chr', 'start', 'pseudobulk_hours']]

    df = pd.merge(df, bulk_rt)

    # compute time from scheduled replication column
    df['time_from_scheduled_rt'] = df['pseudobulk_hours'] - (df[argv.frac_rt_col] * 10.0)

    # find the number of clones and S-phase cells used to compute the Twidth curve
    num_cells = len(df.cell_id.unique())
    num_clones = len(df.assigned_clone_id.unique())

    # dataframe storing computed T-width values
    t_width_df = []

    fig, ax = plt.subplots(1, 2, figsize=(10, 5))

    # compute and plot twidth for per_cell==False
    ax[0], Tw = compute_and_plot_twidth(
        df, tfs_col='time_from_scheduled_rt', rs_col=argv.rep_col, alpha=1,
        title='{} scRT heterogeneity'.format(argv.dataset), curve='sigmoid', ax=ax[0]
    )
    t_width_df.append(pd.DataFrame({
        'dataset': [argv.dataset],
        'per_cell': [False], 'T-width': [Tw],
        'num_cells': [num_cells], 'num_clones': [num_clones]
    }))

    # compute and plot twidth for per_cell==True
    ax[1], Tw = compute_and_plot_twidth(
        df, tfs_col='time_from_scheduled_rt', rs_col=argv.rep_col, alpha=0.1, per_cell=True,
        title='{} scRT heterogeneity'.format(argv.dataset), curve='sigmoid', ax=ax[1]
    )
    t_width_df.append(pd.DataFrame({
        'dataset': [argv.dataset],
        'per_cell': [True], 'T-width': [Tw],
        'num_cells': [num_cells], 'num_clones': [num_clones]
    }))

    # save figure of twidth curves
    fig.savefig(argv.output_pdf, bbox_inches='tight', dpi=300)

    # save a table of all the computed T-width values
    t_width_df = pd.concat(t_width_df, ignore_index=True)
    t_width_df.to_csv(argv.output_csv, index=False)


if __name__ == '__main__':
    main()
