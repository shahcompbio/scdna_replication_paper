import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix, plot_cell_cn_profile
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='matrix of cn pseudobulks')
    p.add_argument('dataset')
    p.add_argument('plot1', type=str, help='pseudobulk cn profiles including dataset-level bulk, sorted by clone_id')
    p.add_argument('plot2', type=str, help='pseudobulk cn profiles with only clones, sorted by clustering')

    return p.parse_args()


def plot_sorted(bulk_df, argv):
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    plot_data = plot_clustered_cell_cn_matrix(
        ax, bulk_df, 'state', cluster_field_name='cluster_id', secondary_field_name='clone_id'
    )
    labels = plot_data.columns.get_level_values(0).values
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title('{}: Pseudobulk CN Profiles'.format(argv.dataset))

    fig.savefig(argv.plot1, bbox_inches='tight')


def plot_clustered(bulk_df, argv):
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    plot_data = plot_clustered_cell_cn_matrix(
        ax, bulk_df.query('cluster_id=="clone"'), 'state', cluster_field_name='cluster_id'
    )
    labels = plot_data.columns.get_level_values(0).values
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels)
    ax.set_title('{}: Pseudobulk CN Profiles'.format(argv.dataset))

    fig.savefig(argv.plot2, bbox_inches='tight')


def main():
    argv = get_args()
    # load long-form dataframes from different cell cycle phases
    bulk_cn = pd.read_csv(argv.input, sep='\t')

    # treat each pseudobulk profile as a unique cell_id
    bulk_df = bulk_cn.melt(id_vars=['chr', 'start'], var_name='cell_id', value_name='state')

    # cluster_id corresponds to whether each profile is dataset, sample or clone level
    bulk_df['cluster_id'] = bulk_df['cell_id'].apply(lambda x: x.split('_')[0])
    # name of the dataset, sample or clone
    bulk_df['clone_id'] = bulk_df['cell_id'].apply(lambda x: x.split('_')[1])

    # plot cn profiles with dataset and sample bulks, sorted by clone_id
    plot_sorted(bulk_df, argv)

    # plot cn profiles without dataset and sample bulks, sorted by clustering
    plot_clustered(bulk_df, argv)


if __name__ == '__main__':
    main()
