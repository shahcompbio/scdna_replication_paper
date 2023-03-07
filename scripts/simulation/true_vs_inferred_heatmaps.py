import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from scgenome import cncluster
from matplotlib.patches import Patch
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_rt_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('dataset', type=str, help='dataset name (e.g. D1.0)')
    p.add_argument('true_frac_col', type=str, help='column name for true fraction of replicated bins per cell')
    p.add_argument('model_rep_state', type=str, help='column name for inferred replication state at each bin')
    p.add_argument('true_rep_state', type=str, help='column name for true replication state at each bin')
    p.add_argument('model_cn_state', type=str, help='column name for inferred copy number state at each bin')
    p.add_argument('true_cn_state', type=str, help='column name for true copy number state at each bin')
    p.add_argument('output', help='heatmap of reads per million, true/inferred replication/CN states, and false positive/negative calls')

    return p.parse_args()


# helper functions for plotting heatmaps
def plot_colorbar(ax, color_mat, title=None):
    """ Given an axis and a color matrix, plot a colorbar. """
    ax.imshow(np.array(color_mat)[::-1, np.newaxis], aspect='auto', origin='lower')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


def plot_color_legend(ax, color_map, title=None):
    """ Given an axis and a color map, plot a legend. """
    legend_elements = []
    for name, color in color_map.items():
        legend_elements.append(Patch(facecolor=color, label=name))
    ax.legend(handles=legend_elements, loc='center left', title=title)
    ax.grid(False)
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])


def make_color_mat_float(values, palette_color):
    """
    Make a color_mat for a 0-1 float array `values` and a
    corresponding color pallete.
    """
    pal = plt.get_cmap(palette_color)
    color_mat = []
    for val in values:
        color_mat.append(pal(val))
    color_dict = {0: pal(0.0), 1: pal(1.0)}
    return color_mat, color_dict


def plot_true_vs_inferred_heatmaps(df, argv):
    """ 
    Plot heatmaps of true and inferred replication and copy states, as well as false positive/negative calls and input read depth.
    All rows of each heatmap are sorted by clone id and the true time in S-phase.
    """
    # create mapping of clones to cluster ids
    cluster_col = 'cluster_id'
    clone_col = 'clone_id'
    clone_dict = dict([(y,x+1) for x,y in enumerate(sorted(df[clone_col].unique()))])
    df[cluster_col] = df[clone_col]
    df = df.replace({cluster_col: clone_dict})

    # note which columns to sort by
    secondary_sort_column = argv.true_frac_col
    secondary_sort_label = 'Time in S-phase'

    # create a color map for the replication states and accuracies
    rt_cmap = get_rt_cmap()

    # compute accuracy of inferred rt_state values
    rep_accuracy = 1.0 - (sum(abs(df[argv.true_rep_state] - df[argv.model_rep_state])) / df.shape[0])

    # compute accuracy of inferred cn states
    cn_accuracy = sum(df[argv.true_cn_state] == df[argv.model_cn_state]) / df.shape[0]

    # create a 2x2 grid of colorbars for the replication and cn states
    fig = plt.figure(figsize=(12, 10))
    ax0 = fig.add_axes([0.05, 0.5, 0.47, 0.45]) # top left
    ax1 = fig.add_axes([0.53, 0.5, 0.47, 0.45]) # top right
    ax2 = fig.add_axes([0.05, 0.0, 0.47, 0.45]) # bottom left
    ax3 = fig.add_axes([0.53, 0.0, 0.47, 0.45]) # bottom right

    # top left: true CN state
    plot_data0 = plot_clustered_cell_cn_matrix(ax0, df, argv.true_cn_state, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column)
    ax0.set_title('{}: True CN state'.format(argv.dataset))

    # top right: inferred CN state
    plot_data1 = plot_clustered_cell_cn_matrix(ax1, df, argv.model_cn_state, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column)
    ax1.set_title('{}: Inferred CN state (accuracy={})'.format(argv.dataset, round(cn_accuracy, 3)))

    # bottom left: true replication state
    plot_data2 = plot_clustered_cell_cn_matrix(ax2, df, argv.true_rep_state, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column, cmap=rt_cmap)
    ax2.set_title('{}: True replication state'.format(argv.dataset))

    # bottom right: inferred replication state
    plot_data3 = plot_clustered_cell_cn_matrix(ax3, df, argv.model_rep_state, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column, cmap=rt_cmap)
    ax3.set_title('{}: Inferred replication state (accuracy={})'.format(argv.dataset, round(rep_accuracy, 3)))

    # hide the y-ticks and labels for all heatmaps
    for ax in [ax0, ax1, ax2, ax3]:
        ax.set_yticks([])
        ax.set_yticklabels([])
    
    if len(clone_dict) > 1:
        # annotate the clones for G1-phase cells
        cell_ids = plot_data0.columns.get_level_values(0).values
        cluster_ids0 = plot_data0.columns.get_level_values(1).values
        color_mat0, color_map0 = cncluster.get_cluster_colors(cluster_ids0, return_map=True)

        # get list of color pigments in the same order as clone_dict
        colors_used0 = []
        for c in color_mat0:
            if c not in colors_used0:
                colors_used0.append(c)

        # match clone IDs to color pigments
        clones_to_colors0 = {}
        for i, key in enumerate(clone_dict.keys()):
            clones_to_colors0[key] = colors_used0[i]

        # get array of secondary_sort_column values that that match the cell_id order
        condensed_cn = df[['cell_id', secondary_sort_column]].drop_duplicates()
        secondary_array = []
        for cell in cell_ids:
            s = condensed_cn[condensed_cn['cell_id'] == cell][secondary_sort_column].values[0]
            secondary_array.append(s)

        # make color mat according to secondary array
        secondary_color_mat, secondary_to_colors = make_color_mat_float(secondary_array, 'Blues')

        # create color bar that shows clone id to the left of the heatmap in the top left corner
        # color bar should be 0.01 wide and 0.45 tall
        ax = fig.add_axes([0.03, 0.5, 0.01, 0.45])
        plot_colorbar(ax, color_mat0)

        # create color bar that shows secondary sort value just to the right of the clone id color bar
        # color bar should be 0.01 wide and 0.45 tall
        ax = fig.add_axes([0.04, 0.5, 0.01, 0.45])
        plot_colorbar(ax, secondary_color_mat)

        # add the same color bars next to the heatmap on the bottom left corner
        ax = fig.add_axes([0.03, 0.0, 0.01, 0.45])
        plot_colorbar(ax, color_mat0)
        ax = fig.add_axes([0.04, 0.0, 0.01, 0.45])
        plot_colorbar(ax, secondary_color_mat)

    # save figure
    fig.savefig(argv.output, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype(str)
    df.chr = df.chr.astype('category')

    # 1 is false positive, 0 is accurate, -1 is false negative
    df['rt_state_diff'] = df[argv.model_rep_state] - df[argv.true_rep_state]

    # create separate plots
    plot_true_vs_inferred_heatmaps(df, argv)


if __name__ == '__main__':
    main()
