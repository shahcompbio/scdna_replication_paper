import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from matplotlib import colors as mcolors
from matplotlib.patches import Patch
from scgenome import cncluster
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_rt_cmap, get_clone_cmap

def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='long-form dataframe of S-phase cells with pyro model results')
    p.add_argument('cn_g', help='long-form dataframe of G1/2-phase cells')
    p.add_argument('dataset')
    p.add_argument('plot1', help='heatmaps of all S-phase cells sorted the same')
    p.add_argument('plot2', help='heatmaps of G1- vs S-phase hmmcopy states')
    p.add_argument('plot3', help='heatmaps of G1- vs S-phase reads per million')

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


def plot_model_results(cn_s, cn_g, argv):
    rt_cmap = get_rt_cmap()
    clone_cmap = get_clone_cmap()

    # create mapping of clone IDs
    cluster_col = 'cluster_id'
    clone_dict = dict([(y,x+1) for x,y in enumerate(sorted(cn_g['clone_id'].unique()))])
    cn_g[cluster_col] = cn_g['clone_id']
    cn_g = cn_g.replace({cluster_col: clone_dict})
    cn_s[cluster_col] = cn_s['clone_id']
    cn_s = cn_s.replace({cluster_col: clone_dict})

    # plot the heatmaps
    fig = plt.figure(figsize=(28,14))

    # plot the S-phase cells in the top row
    # top left corner is the rpm
    ax0 = fig.add_axes([0.05,0.5,0.23,0.45])
    plot_data0 = plot_clustered_cell_cn_matrix(
        ax0, cn_s, 'rpm', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep', max_cn=None, raw=True, cmap='viridis'
    )
    ax0.set_title('{} S-phase cells\nReads per million'.format(argv.dataset))

    # top mid-left is the hmmcopy states
    ax1 = fig.add_axes([0.29,0.5,0.23,0.45])
    plot_data1 = plot_clustered_cell_cn_matrix(
        ax1, cn_s, 'state', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep'
    )
    ax1.set_title('{} S-phase cells\nHMMcopy states'.format(argv.dataset))

    # top mid-right is the model cn states
    ax2 = fig.add_axes([0.53,0.5,0.23,0.45])
    plot_data2 = plot_clustered_cell_cn_matrix(
        ax2, cn_s, 'model_cn_state', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep'
    )
    ax2.set_title('{} S-phase cells\nPERT CN states'.format(argv.dataset))

    # top right corner is the replication states
    ax3 = fig.add_axes([0.77,0.5,0.23,0.45])
    plot_data3 = plot_clustered_cell_cn_matrix(
        ax3, cn_s, 'model_rep_state', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep', cmap=rt_cmap
    )
    ax3.set_title('{} S-phase cells\nPERT replication states'.format(argv.dataset))

    # plot the G1/2-phase cells in the bottom row
    # bottom left corner is the rpm
    ax4 = fig.add_axes([0.05,0.0,0.23,0.45])
    plot_data4 = plot_clustered_cell_cn_matrix(
        ax4, cn_g, 'rpm', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep', max_cn=None, raw=True, cmap='viridis'
    )
    ax4.set_title('{} G1/2-phase cells\nReads per million'.format(argv.dataset))

    # bottom mid-left is the hmmcopy states
    ax5 = fig.add_axes([0.29,0.0,0.23,0.45])
    plot_data5 = plot_clustered_cell_cn_matrix(
        ax5, cn_g, 'state', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep'
    )
    ax5.set_title('{} G1/2-phase cells\nHMMcopy states'.format(argv.dataset))

    # bottom mid-right is the model cn states
    ax6 = fig.add_axes([0.53,0.0,0.23,0.45])
    plot_data6 = plot_clustered_cell_cn_matrix(
        ax6, cn_g, 'model_cn_state', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep'
    )
    ax6.set_title('{} G1/2-phase cells\nPERT CN states'.format(argv.dataset))

    # bottom right corner is the replication states
    ax7 = fig.add_axes([0.77,0.0,0.23,0.45])
    plot_data7 = plot_clustered_cell_cn_matrix(
        ax7, cn_g, 'model_rep_state', cluster_field_name=cluster_col, secondary_field_name='cell_frac_rep', cmap=rt_cmap
    )
    ax7.set_title('{} G1/2-phase cells\nPERT replication states'.format(argv.dataset))

    # turn off the y-axis ticks in all subplots
    for ax in [ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7]:
        ax.set_yticks([])
        ax.set_ylabel('')

    # add the colorbars for clone_id and cell_frac_rep
    if len(clone_dict) > 1:
        cell_ids = plot_data0.columns.get_level_values(0).values
        cluster_ids0 = plot_data0.columns.get_level_values(1).values
        # use mcolors to change every element in the dict to rgba
        for key in clone_cmap.keys():
            clone_cmap[key] = mcolors.to_rgba(clone_cmap[key])
        color_mat0, color_map0 = cncluster.get_cluster_colors(cluster_ids0, color_map=clone_cmap, return_map=True)

        # get array of 'cell_frac_rep' values that that match the cell_id order
        condensed_cn = cn_s[['cell_id', 'cell_frac_rep']].drop_duplicates()
        secondary_array = []
        for cell in cell_ids:
            s = condensed_cn[condensed_cn['cell_id'] == cell]['cell_frac_rep'].values[0]
            secondary_array.append(s)

        # make color mat according to secondary array
        secondary_color_mat, secondary_to_colors = make_color_mat_float(secondary_array, 'Blues')

        # create color bar that shows clone id for each row in heatmap
        ax = fig.add_axes([0.03,0.5,0.01,0.45])
        plot_colorbar(ax, color_mat0)

        # create color bar that shows secondary sort value for each row in heatmap
        ax = fig.add_axes([0.04,0.5,0.01,0.45])
        plot_colorbar(ax, secondary_color_mat)

        # repeat for the G1/2-phase cells in the bottom row
        cell_ids = plot_data4.columns.get_level_values(0).values
        cluster_ids4 = plot_data4.columns.get_level_values(1).values
        color_mat4, color_map4 = cncluster.get_cluster_colors(cluster_ids4, color_map=clone_cmap, return_map=True)

        # get array of 'cell_frac_rep' values that that match the cell_id order
        condensed_cn = cn_g[['cell_id', 'cell_frac_rep']].drop_duplicates()
        secondary_array = []
        for cell in cell_ids:
            s = condensed_cn[condensed_cn['cell_id'] == cell]['cell_frac_rep'].values[0]
            secondary_array.append(s)
        
        # make color mat according to secondary array
        secondary_color_mat, secondary_to_colors = make_color_mat_float(secondary_array, 'Blues')

        # create color bar that shows clone id for each row in heatmap
        ax = fig.add_axes([0.03,0.0,0.01,0.45])
        plot_colorbar(ax, color_mat4)

        # create color bar that shows secondary sort value for each row in heatmap
        ax = fig.add_axes([0.04,0.0,0.01,0.45])
        plot_colorbar(ax, secondary_color_mat)

    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def plot_hmmcopy(cn_s, cn_g1, argv):
    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    plot_clustered_cell_cn_matrix(ax[0], cn_g1, 'state', cluster_field_name='clone_id')
    plot_clustered_cell_cn_matrix(ax[1], cn_s, 'state', cluster_field_name='clone_id')

    ax[0].set_title('{}\nG1-phase HMMCopy states'.format(argv.dataset))
    ax[1].set_title('{}\nS-phase HMMCopy states'.format(argv.dataset))
    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)


def plot_rpm(cn_s, cn_g1, argv):
    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    plot_clustered_cell_cn_matrix(ax[0], cn_g1, 'rpm', max_cn=None, raw=True, cmap='viridis', cluster_field_name='clone_id')
    plot_clustered_cell_cn_matrix(ax[1], cn_s, 'rpm', max_cn=None, raw=True, cmap='viridis', cluster_field_name='clone_id')

    ax[0].set_title('{}\nG1-phase reads per million'.format(argv.dataset))
    ax[1].set_title('{}\nS-phase reads per million'.format(argv.dataset))
    fig.savefig(argv.plot3, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g = pd.read_csv(argv.cn_g, sep='\t')

    # show rpm, hmmcopy, inferred cn, inferred rep heatmaps for S-phase cells and G1/2-phase cells
    # where all the rows are sorted the same in all four heatmaps
    plot_model_results(cn_s, cn_g, argv)

    # show hmmcopy state heatmaps for both S-phase and G1-phase cells
    plot_hmmcopy(cn_s, cn_g, argv)

    # show reads per million heatmaps for both S-phase and G1-phase cells
    plot_rpm(cn_s, cn_g, argv)



if __name__=='__main__':
    main()
