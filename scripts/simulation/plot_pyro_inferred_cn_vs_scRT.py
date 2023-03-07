from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
from scgenome.cnplot import plot_clustered_cell_cn_matrix
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_rt_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with inferred cn and replication states for each bin')
    p.add_argument('rep_col', help='column name for inferred replication states')
    p.add_argument('cn_col', help='column name for inferred copy number')
    p.add_argument('frac_rt_col', help='column name for time within S-phase (i.e. fraction replicated)')
    p.add_argument('output_rt_state', help='heatmaps of inferred cn and scRT states')

    return p.parse_args()


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state'):
    ''' Compute the fraction of replicated bins for all cells in `cn` '''
    cn['extreme_cell_frac'] = False
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn[rep_state_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac
    return cn


def plot_cn_and_rep_states(df, argv):

    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    rt_cmap = get_rt_cmap()
    plot_clustered_cell_cn_matrix(ax[0], df, argv.cn_col, cluster_field_name='clone_id', secondary_field_name=argv.frac_rt_col)
    plot_clustered_cell_cn_matrix(ax[1], df, argv.rep_col, cluster_field_name='clone_id', secondary_field_name=argv.frac_rt_col, cmap=rt_cmap)

    ax[0].set_title('Inferred CN states')
    ax[1].set_title('Inferred replication states')

    fig.savefig(argv.output_rt_state, bbox_inches='tight')


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # compute fraction of replicated bins per cells 
    if argv.frac_rt_col not in df.columns:
        df = compute_cell_frac(df, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.rep_col)

    # create separate plots
    plot_cn_and_rep_states(df, argv)


if __name__ == '__main__':
    main()
