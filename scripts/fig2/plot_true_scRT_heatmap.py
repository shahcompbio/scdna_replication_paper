from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from scgenome import cncluster
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true scRT data')
    p.add_argument('dataset')
    p.add_argument('output', help='heatmap of read count and true scRT states')

    return p.parse_args()


# helper functions for plotting heatmaps
def plot_colorbar(ax, color_mat, title=None):
    ax.imshow(np.array(color_mat)[::-1, np.newaxis], aspect='auto', origin='lower')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


def plot_color_legend(ax, color_map, title=None):
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


def get_rt_cmap():
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list)


def plot_true_rt_state(df, argv):
    df = df.copy()

    # create mapping of clones
    cluster_col = 'cluster_id'
    clone_col = 'clone_id'
    clone_dict = dict([(y,x+1) for x,y in enumerate(sorted(df[clone_col].unique()))])
    df[cluster_col] = df[clone_col]
    df = df.replace({cluster_col: clone_dict})

    secondary_sort_column = 'true_t'
    secondary_sort_label = 'Time in S-phase'

    fig = plt.figure(figsize=(14, 7))
    
    ax0 = fig.add_axes([0.12,0.0,0.38,1.])
    rt_cmap = get_rt_cmap()
    plot_data0 = plot_clustered_cell_cn_matrix(ax0, df, 'true_rep', cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column, cmap=rt_cmap)
    ax0.set_title('{}: True rep state'.format(argv.dataset))

    ax1 = fig.add_axes([0.62,0.0,0.38,1.])
    plot_data1 = plot_clustered_cell_cn_matrix(ax1, df, 'true_reads_norm', cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column, cmap='viridis', max_cn=None)
    ax1.set_title('{}: Read count'.format(argv.dataset))
    
    
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

        # create color bar that shows clone id for each row in heatmap
        ax = fig.add_axes([0.09,0.0,0.03,1.])
        plot_colorbar(ax, color_mat0)

        # create color bar that shows secondary sort value for each row in heatmap
        ax = fig.add_axes([0.06,0.0,0.03,1.])
        plot_colorbar(ax, secondary_color_mat)

        # create legend to match colors to clone ids
        ax = fig.add_axes([0.0,0.75,0.04,0.25])
        plot_color_legend(ax, clones_to_colors0, title='Clone ID')

        # create legend to match colors to secondary sort values
        ax = fig.add_axes([0.0,0.5,0.04,0.25])
        plot_color_legend(ax, secondary_to_colors, title=secondary_sort_label)

        # annotate the clones for S-phase cells.. using the same colors as G1 clones
        cluster_ids1 = plot_data1.columns.get_level_values(1).values
        color_mat1 = cncluster.get_cluster_colors(cluster_ids1, color_map=color_map0)

        # match clone IDs to color pigments
        clones_to_colors1 = {}
        for i, key in enumerate(clone_dict.keys()):
            clones_to_colors1[key] = colors_used0[i]

        # create color bar that shows clone id for each row in heatmap
        ax = fig.add_axes([0.59,0.0,0.03,1.])
        plot_colorbar(ax, color_mat1)

        # create color bar that shows secondary sort value for each row in heatmap
        ax = fig.add_axes([0.56,0.0,0.03,1.])
        plot_colorbar(ax, secondary_color_mat)

        # create legend to match colors to clone ids
        ax = fig.add_axes([0.5,0.75,0.04,0.25])
        plot_color_legend(ax, clones_to_colors1, title='Clone ID')

         # create legend to match colors to secondary sort values
        ax = fig.add_axes([0.5,0.5,0.04,0.25])
        plot_color_legend(ax, secondary_to_colors, title=secondary_sort_label)

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
