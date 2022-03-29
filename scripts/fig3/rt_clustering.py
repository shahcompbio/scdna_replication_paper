from argparse import ArgumentParser
import umap
import hdbscan
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
    p.add_argument('value_col', help='column containing rt states')
    p.add_argument('sort_col', help='column containing time within S-phase')
    p.add_argument('dataset')
    p.add_argument('output_umap', help='heatmap of inferred scRT states')
    p.add_argument('output_kde', help='kdeplot of pct replicated split by cluster')
    p.add_argument('output_heatmap', help='heatmap of inferred scRT states split by cluster')
    p.add_argument('output_df', help='cn_s input df with scRT clusters added')

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


def get_rt_state_clusters(rt):
    # generate UMAP embedding
    reducer = umap.UMAP()
    embedding = reducer.fit_transform(rt.values)    

    # run hdbscan to find clusters
    clusterer = hdbscan.HDBSCAN()
    clusterer.fit(embedding)

    return embedding, clusterer


def plot_umaps(rt, embedding, clusterer, argv):
    fig, ax = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=rt.reset_index()[argv.sort_col].values, ax=ax[0])
    sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=rt.reset_index()['clone_id'].values, ax=ax[1])
    sns.scatterplot(x=embedding[:, 0], y=embedding[:, 1], hue=clusterer.labels_.astype(str), ax=ax[2])

    ax[0].set_title('UMAP projection of {}'.format(argv.dataset))
    ax[1].set_title('UMAP projection of {}'.format(argv.dataset))
    ax[2].set_title('UMAP projection of {}'.format(argv.dataset))
    ax[0].legend(title='% Rep')
    ax[1].legend(title='Clone ID')
    ax[2].legend(title='RT cluster')

    fig.savefig(argv.output_umap, bbox_inches='tight')


def plot_kde(summary, argv):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    sns.kdeplot(data=summary, x=argv.sort_col, hue='hdbscan_cluster', ax=ax)
    ax.set_xlabel('% Replicated')
    ax.set_title(argv.dataset)
    fig.savefig(argv.output_kde, bbox_inches='tight')


def plot_heatmaps(df, clusterer, argv):
    N = clusterer.labels_.max() - clusterer.labels_.min() + 1
    print('N', N)
    fig, ax = plt.subplots(N, 1, figsize=(10, 5*N), tight_layout=True)
    ax = ax.flatten()
    rt_cmap = get_rt_cmap()

    df['dummy_id'] = 0
    i = 0
    for rt_cluster, chunk in df.groupby('hdbscan_cluster'):
        print('i', i)
        print('rt_cluster', rt_cluster)
        _ = plot_clustered_cell_cn_matrix(ax[i], chunk, argv.value_col, cluster_field_name='dummy_id', secondary_field_name=argv.sort_col, cmap=rt_cmap)
        ax[i].set_title('{}, RT cluster {}'.format(argv.dataset, rt_cluster))
        i += 1

    fig.savefig(argv.output_heatmap, bbox_inches='tight')



def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype(str)
    df.chr = df.chr.astype('category')

    # pivot into RT state matrices
    rt = df.pivot(index=['cell_id', argv.sort_col, 'clone_id'], columns=['chr', 'start'], values=argv.value_col)

    # drop loci that might only be present in a handful of cells
    rt = rt.dropna(axis=1)

    embedding, clusterer = get_rt_state_clusters(rt)
    plot_umaps(rt, embedding, clusterer, argv)

    # merge into one df
    summary = pd.DataFrame({'cell_id': rt.reset_index()['cell_id'].values, 
                           'clone_id': rt.reset_index()['clone_id'].values,
                           argv.sort_col: rt.reset_index()[argv.sort_col].values,
                           'hdbscan_cluster': clusterer.labels_.astype(str)
                           })
    df = pd.merge(df, summary)

    # show kdeplot of % replicated to further illustrate that most clusters
    # are from different times within S-phase
    plot_kde(summary, argv)

    # plot clustered heatmap for each state
    plot_heatmaps(df, clusterer, argv)

    df.to_csv(argv.output_df, sep='\t', index=False)


if __name__ == '__main__':
    main()
