from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scdna_replication_tools.calculate_twidth import compute_and_plot_twidth


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('bulk_rt', help='pseduobulk RT information')
    p.add_argument('dataset')
    p.add_argument('output_pdf', help='T-width curve using inferred scRT')

    return p.parse_args()


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    bulk_rt = pd.read_csv(argv.bulk_rt, sep='\t')
    bulk_rt = bulk_rt[['chr', 'start', 'pseduobulk_rt_value', 'pseduobulk_hours']]

    df = pd.merge(df, bulk_rt)

    # compute time from scheduled replication column
    df['time_from_scheduled_rt'] = df['pseduobulk_hours'] - (df['frac_rt'] * 10.0)

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    ax, t_width = compute_and_plot_twidth(df, tfs_col='time_from_scheduled_rt', rs_col='rt_state', alpha=1,
                                          title='{} scRT heterogeneity'.format(argv.dataset), curve='sigmoid', ax=ax)

    fig.savefig(argv.output_pdf, bbox_inches='tight')


if __name__ == '__main__':
    main()
