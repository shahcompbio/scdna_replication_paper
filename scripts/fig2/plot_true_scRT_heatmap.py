from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from matplotlib.colors import ListedColormap


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true scRT data')
    p.add_argument('dataset')
    p.add_argument('output', help='heatmap of read count and true scRT states')

    return p.parse_args()


def get_rt_cmap():
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list)


def plot_true_rt_state(df, argv):

    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    rt_cmap = get_rt_cmap()
    plot_clustered_cell_cn_matrix(ax[0], df, 'true_rt_state', cluster_field_name='clone_id', secondary_field_name='true_frac_rt', cmap=rt_cmap)
    ax[0].set_title('{}: True scRT'.format(argv.dataset))

    plot_clustered_cell_cn_matrix(ax[1], df, 'reads', cluster_field_name='clone_id', secondary_field_name='true_frac_rt', cmap='viridis')
    ax[1].set_title('{}: Read count'.format(argv.dataset))

    fig.savefig(argv.output, bbox_inches='tight')


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype(str)
    df.chr = df.chr.astype('category')

    # create separate plots
    plot_true_rt_state(df, argv)


if __name__ == '__main__':
    main()
