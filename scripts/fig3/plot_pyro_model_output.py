import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix, plot_cell_cn_profile
from matplotlib.colors import ListedColormap
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='long-form dataframe of S-phase cells with pyro model results')
    p.add_argument('cn_g', help='long-form dataframe of G1/2-phase cells')
    p.add_argument('dataset')
    p.add_argument('plot1', help='heatmaps of all S-phase cells sorted the same')
    p.add_argument('plot2', help='heatmaps of G1- vs S-phase hmmcopy states')
    p.add_argument('plot3', help='heatmaps of G1- vs S-phase reads per million')

    return p.parse_args()


def get_rt_cmap():
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list)


def plot_model_results(cn_s, argv):
    rt_cmap = get_rt_cmap()

    fig, ax = plt.subplots(1, 4, figsize=(28, 7), tight_layout=True)
    ax = ax.flatten()

    plot_clustered_cell_cn_matrix(ax[0], cn_s, 'rpm', max_cn=None, raw=True, cmap='viridis', cluster_field_name='clone_id', secondary_field_name='model_s_time')
    plot_clustered_cell_cn_matrix(ax[1], cn_s, 'state', cluster_field_name='clone_id', secondary_field_name='model_s_time')
    plot_clustered_cell_cn_matrix(ax[2], cn_s, 'model_cn_state', cluster_field_name='clone_id', secondary_field_name='model_s_time')
    plot_clustered_cell_cn_matrix(ax[3], cn_s, 'model_rep_state', cluster_field_name='clone_id', secondary_field_name='model_s_time', cmap=rt_cmap)

    ax[0].set_title('{}\nReads per million'.format(argv.dataset))
    ax[1].set_title('{}\nHMMCopy states'.format(argv.dataset))
    ax[2].set_title('{}\nInferred CN states'.format(argv.dataset))
    ax[3].set_title('{}\nInferred replication states'.format(argv.dataset))
    fig.savefig(argv.plot1, bbox_inches='tight')


def plot_hmmcopy(cn_s, cn_g1, argv):
    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    plot_clustered_cell_cn_matrix(ax[0], cn_g1, 'state', cluster_field_name='clone_id')
    plot_clustered_cell_cn_matrix(ax[1], cn_s, 'state', cluster_field_name='clone_id')

    ax[0].set_title('{}\nG1-phase HMMCopy states'.format(argv.dataset))
    ax[1].set_title('{}\nS-phase HMMCopy states'.format(argv.dataset))
    fig.savefig(argv.plot2, bbox_inches='tight')


def plot_rpm(cn_s, cn_g1, argv):
    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    plot_clustered_cell_cn_matrix(ax[0], cn_g1, 'rpm', max_cn=None, raw=True, cmap='viridis', cluster_field_name='clone_id')
    plot_clustered_cell_cn_matrix(ax[1], cn_s, 'rpm', max_cn=None, raw=True, cmap='viridis', cluster_field_name='clone_id')

    ax[0].set_title('{}\nG1-phase reads per million'.format(argv.dataset))
    ax[1].set_title('{}\nS-phase reads per million'.format(argv.dataset))
    fig.savefig(argv.plot3, bbox_inches='tight')


def main():
    argv = get_args()

    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g = pd.read_csv(argv.cn_g, sep='\t')

    # show rpm, hmmcopy, inferred cn, inferred rep heatmaps for S-phase cells
    # where all the rows are sorted the same in all four heatmaps
    plot_model_results(cn_s, argv)

    # show hmmcopy state heatmaps for both S-phase and G1-phase cells
    plot_hmmcopy(cn_s, cn_g, argv)

    # show reads per million heatmaps for both S-phase and G1-phase cells
    plot_rpm(cn_s, cn_g, argv)



if __name__=='__main__':
    main()
