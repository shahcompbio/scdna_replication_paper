import matplotlib
matplotlib.use('Agg')
import sys
import scgenome.cnplot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from scgenome import cncluster
from matplotlib.patches import Patch


def get_args():
    p = ArgumentParser()

    p.add_argument('cn', help='cn input tsv file for cells of all input phases')
    p.add_argument('value_col')
    p.add_argument('dataset')
    p.add_argument('output_png', help='plot containinig copy number heatmaps')

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


def plot_cn_heatmap(cn_g, cn_s, figsize=(14,7), dataset=None, clone_col='library_id', value_col='state'):
    cn_g = cn_g.copy()
    cn_s = cn_s.copy()

    # create mapping of cell clones
    cluster_col = 'cluster_id'
    clone_dict = dict([(y,x+1) for x,y in enumerate(sorted(cn_g[clone_col].unique()))])
    cn_g[cluster_col] = cn_g[clone_col]
    cn_g = cn_g.replace({cluster_col: clone_dict})
    cn_s[cluster_col] = cn_s[clone_col]
    cn_s = cn_s.replace({cluster_col: clone_dict})

    # plot the heatmaps
    fig = plt.figure(figsize=figsize)
    ax_g1 = fig.add_axes([0.1,0.0,0.425,1.])
    plot_data_g1 = scgenome.cnplot.plot_clustered_cell_cn_matrix(
        ax_g1, cn_g, value_col, cluster_field_name=cluster_col
    )

    ax_s = fig.add_axes([0.575,0.0,0.425,1.])
    plot_data_s = scgenome.cnplot.plot_clustered_cell_cn_matrix(
        ax_s, cn_s, value_col, cluster_field_name=cluster_col
    )

    # turn off the y-axis ticks
    for ax in [ax_g1, ax_s]:
        ax.set_yticks([])
        ax.set_ylabel('')
    
    if dataset:
        ax_g1.set_title('{}: G1/2-phase high quality PERT input'.format(dataset))
        ax_s.set_title('{}: S-phase & low quality PERT input'.format(dataset))
    
    if len(clone_dict) > 1:
        # annotate the clones for G1-phase cells
        cluster_ids_g1 = plot_data_g1.columns.get_level_values(1).values
        color_mat_g1, color_map_g1 = cncluster.get_cluster_colors(cluster_ids_g1, return_map=True)

        # get list of color pigments in the same order as clone_dict
        colors_used_g1 = []
        for c in color_mat_g1:
            if c not in colors_used_g1:
                colors_used_g1.append(c)

        # match clone IDs to color pigments
        clones_to_colors_g1 = {}
        for i, key in enumerate(clone_dict.keys()):
            clones_to_colors_g1[key] = colors_used_g1[i]

        # create color bar that shows clone id for each row in heatmap
        ax = fig.add_axes([0.07,0.0,0.03,1.])
        plot_colorbar(ax, color_mat_g1)

        # create legend to match colors to clone ids
        ax = fig.add_axes([0.0,0.75,0.02,0.25])
        plot_color_legend(ax, clones_to_colors_g1, title='Clone')

        # # annotate the clones for S-phase cells.. using the same colors as G1 clones
        # cluster_ids_s = plot_data_s.columns.get_level_values(1).values
        # color_mat_s = cncluster.get_cluster_colors(cluster_ids_s, color_map=color_map_g1)

        # # match clone IDs to color pigments
        # clones_to_colors_s = {}
        # for i, key in enumerate(clone_dict.keys()):
        #     clones_to_colors_s[key] = colors_used_g1[i]

        # # create color bar that shows clone id for each row in heatmap
        # ax = fig.add_axes([0.545,0.0,0.03,1.])
        # plot_colorbar(ax, color_mat_s)


    return fig



def main():
    argv = get_args()

    cn = pd.read_csv(argv.cn, sep='\t')

    # convert the 'chr' column to a string and then categorical
    cn.chr = cn.chr.astype(str)
    cn.chr = cn.chr.astype('category')

    # split into G1/2 and S phases based on the 'clonealign_clone_id' column
    cn_g = cn[cn['clone_id']!='not_in_tree']
    cn_s = cn[cn['clone_id']=='not_in_tree']

    fig = plot_cn_heatmap(cn_g, cn_s, dataset=argv.dataset, value_col=argv.value_col, clone_col='clone_id')

    fig.savefig(argv.output_png, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
