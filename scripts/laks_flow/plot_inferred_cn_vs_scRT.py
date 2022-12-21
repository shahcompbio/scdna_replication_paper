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

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with inferred cn and replication states for each bin')
    p.add_argument('rep_col', help='column name for inferred replication states')
    p.add_argument('cn_col', help='column name for inferred copy number')
    p.add_argument('frac_rt_col', help='column name for time within S-phase (i.e. fraction replicated)')
    p.add_argument('output_rt_state', help='heatmap comparing true and inferred rt_state values')
    p.add_argument('output_frac_rt', help='plot for frac_rt distribution across all cells')

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


def get_rt_cmap():
    """ Return a colormap for replication states. """
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list), rt_colors


def plot_cn_and_rep_states(df, argv):
    """ 
    Plot 3 heatmaps for S-phase cells: hmmcopy state, inferred CN state, inferred rep state
    All rows should be sorted first by clone_id and then by frac_rt. The colorbars should be on the far left.
    """
    # create mapping of clones to cluster ids
    cluster_col = 'cluster_id'
    clone_col = 'clone_id'
    clone_dict = dict([(y,x+1) for x,y in enumerate(sorted(df[clone_col].unique()))])
    df[cluster_col] = df[clone_col]
    df = df.replace({cluster_col: clone_dict})

    # note which columns to sort by
    secondary_sort_column = argv.frac_rt_col
    secondary_sort_label = 'Time in S-phase'

    # create a color map for the replication states and accuracies
    rt_cmap, rt_color_dict = get_rt_cmap()
    
    # add 3 subplots, leaving space for the colorbars on the far left
    # the heatmaps should be arranged in 1 row and 3 columns, having a height of 1 and a width of 0.29
    fig = plt.figure(figsize=(21, 7))
    ax0 = fig.add_axes([0.1, 0.0, 0.29, 1.0]) # left
    ax1 = fig.add_axes([0.4, 0.0, 0.29, 1.0]) # middle
    ax2 = fig.add_axes([0.7, 0.0, 0.29, 1.0]) # right

    # left: hmmcopy states
    plot_data0 = plot_clustered_cell_cn_matrix(ax0, df, 'state', cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column)
    ax0.set_title('HMMcopy states')

    # middle: inferred CN states
    plot_data1 = plot_clustered_cell_cn_matrix(ax1, df, argv.cn_col, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column)
    ax1.set_title('Inferred CN states')

    # right: inferred replication states
    plot_data2 = plot_clustered_cell_cn_matrix(ax2, df, argv.rep_col, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column, cmap=rt_cmap)
    ax2.set_title('Inferred replication states')

    # hide the y-ticks and labels for all heatmaps
    for ax in [ax0, ax1, ax2]:
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
        # color bar should be 0.01 wide and 1.0 tall
        ax = fig.add_axes([0.08, 0.0, 0.01, 1.0])
        plot_colorbar(ax, color_mat0)

        # create color bar that shows secondary sort value just to the right of the clone id color bar
        # color bar should be 0.01 wide and 1.0 tall
        ax = fig.add_axes([0.09, 0.0, 0.01, 1.0])
        plot_colorbar(ax, secondary_color_mat)

        # create legend to match colors to clone ids
        ax = fig.add_axes([0.0, 0.75, 0.04, 0.25])
        plot_color_legend(ax, clones_to_colors0, title='Clone ID')

        # create legend to match colors to secondary sort values
        ax = fig.add_axes([0.0,0.5,0.04,0.25])
        plot_color_legend(ax, secondary_to_colors, title=secondary_sort_label)

        # create a legend for the true and inferred replication states
        ax = fig.add_axes([0.0, 0.25, 0.04, 0.25])
        plot_color_legend(ax, rt_color_dict, title='Rep state')

    # save figure
    fig.savefig(argv.output_rt_state, bbox_inches='tight', dpi=300)


def plot_frac_rt_distributions(df, argv):
    df_frac = df[['cell_id', argv.frac_rt_col, 'library_id', 'clone_id']].drop_duplicates().reset_index(drop=True)
    
    fig, ax = plt.subplots(1,3, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    # violinplot
    sns.histplot(data=df_frac, x=argv.frac_rt_col, ax=ax[0])
    sns.histplot(data=df_frac, x=argv.frac_rt_col, hue='library_id', multiple='stack', ax=ax[1])
    sns.histplot(data=df_frac, x=argv.frac_rt_col, hue='clone_id', multiple='stack', ax=ax[2])

    for i in range(3):
        ax[i].set_xlabel('Time in S-phase')
        ax[i].set_title('Distribution of cells\nwithin S-phase')

    fig.savefig(argv.output_frac_rt, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # create separate plots
    plot_cn_and_rep_states(df, argv)
    plot_frac_rt_distributions(df, argv)


if __name__ == '__main__':
    main()
