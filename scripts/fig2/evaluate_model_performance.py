from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from matplotlib.colors import ListedColormap


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('output_rt_state', help='heatmap comparing true and inferred rt_state values')
    p.add_argument('output_rt_accuracy', help='heatmap of false pos/neg and changepoint segments')
    p.add_argument('output_frac_rt', help='Compare distributions of frac_rt between true and inferred')

    return p.parse_args()


def get_rt_cmap():
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list)


def get_chng_cmap(max_int):
    changepoint_colors = {0:'#CCCCCC'}
    color_list = []
    for i in [0]:
        color_list.append(changepoint_colors[i])
    for i in range(1, max_int+1):
        color_list.append('C{}'.format(i))
    return ListedColormap(color_list)


def get_acc_cmap():
    acc_colors = {0:'#CCCCCC', -1: '#532A44', 1: '#00685E'}
    color_list = []
    for i in [-1, 0, 1]:
        color_list.append(acc_colors[i])
    return ListedColormap(color_list)


def plot_true_vs_inferred_rt_state(df, argv):
    # compute accuracy of inferred rt_state values
    accuracy = 1.0 - (sum(abs(df['true_rt_state'] - df['rt_state'])) / df.shape[0])

    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    rt_cmap = get_rt_cmap()
    plot_clustered_cell_cn_matrix(ax[0], df, 'true_rt_state', secondary_field_name='true_frac_rt', cmap=rt_cmap)
    plot_clustered_cell_cn_matrix(ax[1], df, 'rt_state', secondary_field_name='true_frac_rt', cmap=rt_cmap)

    ax[0].set_title('True scRT')
    ax[1].set_title('Inferred scRT\nAccuracy: {}'.format(round(accuracy, 3)))

    fig.savefig(argv.output_rt_state, bbox_inches='tight')


def plot_rt_accuracy(df, argv):
    accuracy = 1.0 - (sum(abs(df['true_rt_state'] - df['rt_state'])) / df.shape[0])

    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    acc_cmap = get_acc_cmap()
    chng_cmap = get_chng_cmap(int(max(df['changepoint_segments'])))
    plot_clustered_cell_cn_matrix(ax[0], df, 'rt_state_diff', secondary_field_name='true_frac_rt', cmap=acc_cmap)
    plot_clustered_cell_cn_matrix(ax[1], df, 'changepoint_segments', secondary_field_name='true_frac_rt', cmap=chng_cmap)

    ax[0].set_title('False positives (green) and negatives (purple)\nAccuracy: {}'.format(round(accuracy, 3)))
    ax[1].set_title('Inferred changepoint segments')

    fig.savefig(argv.output_rt_accuracy, bbox_inches='tight')


def plot_frac_rt_distributions(df, argv):
    df_frac = df[['cell_id', 'true_frac_rt', 'frac_rt']].drop_duplicates().reset_index(drop=True)
    df_frac2 = df_frac.rename(columns={'true_frac_rt': 'true', 'frac_rt': 'inferred'})
    df_frac2 = df_frac2.melt(id_vars='cell_id', var_name='source', value_name='frac_rt')

    fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True)
    ax = ax.flatten()

    # violinplot
    sns.violinplot(data=df_frac2, x='source', y='frac_rt', ax=ax[0])

    # swarmplot with lines connecting point pairs
    sns.swarmplot(data=df_frac2, x='source', y='frac_rt', ax=ax[1], size=3)
    # Now connect the dots
    # Find idx0 and idx1 by inspecting the elements return from ax.get_children()
    # ... or find a way to automate it
    idx0 = 0
    idx1 = 1
    locs1 = ax[1].get_children()[idx0].get_offsets()
    locs2 = ax[1].get_children()[idx1].get_offsets()

    # before plotting, we need to sort so that the data points
    # correspond to each other as they did in "set1" and "set2"
    set1 = df_frac['true_frac_rt'].values
    set2 = df_frac['frac_rt'].values
    sort_idxs1 = np.argsort(set1)
    sort_idxs2 = np.argsort(set2)

    # revert "ascending sort" through sort_idxs2.argsort(),
    # and then sort into order corresponding with set1
    locs2_sorted = locs2[sort_idxs2.argsort()][sort_idxs1]

    for i in range(locs1.shape[0]):
        x = [locs1[i, 0], locs2_sorted[i, 0]]
        y = [locs1[i, 1], locs2_sorted[i, 1]]
        ax[1].plot(x, y, color="black", alpha=0.05)

    for i in range(2):
        ax[i].set_ylabel('Fraction replicated')
        ax[i].set_title('Distribution of cells\nwithin S-phase')

    fig.savefig(argv.output_frac_rt, bbox_inches='tight')


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # assign dummy cluster_id column for scgenome heatmaps
    df['cluster_id'] = 'A'

    # 1 is false positive, 0 is accurate, -1 is false negative
    df['rt_state_diff'] = df['rt_state'] - df['true_rt_state']

    # create separate plots
    plot_true_vs_inferred_rt_state(df, argv)
    plot_rt_accuracy(df, argv)
    plot_frac_rt_distributions(df, argv)


if __name__ == '__main__':
    main()
