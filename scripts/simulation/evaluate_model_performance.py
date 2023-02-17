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
    p.add_argument('rep_col', help='column name for inferred replication states')
    p.add_argument('cn_col', help='column name containing copy number (pyro model) or changepoint segment (heuristic model)')
    p.add_argument('frac_rt_col', help='column name for time within S-phase (i.e. fraction replicated)')
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


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state'):
    ''' Compute the fraction of replicated bins for all cells in `cn` '''
    cn['extreme_cell_frac'] = False
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn[rep_state_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac
    return cn


def plot_true_vs_inferred_rt_state(df, argv):
    # compute accuracy of inferred rt_state values
    accuracy = 1.0 - (sum(abs(df['true_rep'] - df[argv.rep_col])) / df.shape[0])

    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    rt_cmap = get_rt_cmap()
    plot_clustered_cell_cn_matrix(ax[0], df, 'true_rep', cluster_field_name='clone_id', secondary_field_name='true_t', cmap=rt_cmap)
    plot_clustered_cell_cn_matrix(ax[1], df, argv.rep_col, cluster_field_name='clone_id', secondary_field_name='true_t', cmap=rt_cmap)

    ax[0].set_title('True scRT')
    ax[1].set_title('Inferred scRT\nAccuracy: {}'.format(round(accuracy, 3)))

    fig.savefig(argv.output_rt_state, bbox_inches='tight')


def plot_rt_accuracy(df, argv):
    accuracy = 1.0 - (sum(abs(df['true_rep'] - df[argv.rep_col])) / df.shape[0])

    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    if argv.cn_col=='observed_cn_state':
        plot_clustered_cell_cn_matrix(ax[0], df, argv.cn_col, cluster_field_name='clone_id', secondary_field_name='true_t')
        ax[0].set_title('True CN state')
    elif 'cn' in argv.cn_col:
        plot_clustered_cell_cn_matrix(ax[0], df, argv.cn_col, cluster_field_name='clone_id', secondary_field_name='true_t')
        ax[0].set_title('Inferred CN state')
    else:
        print('specify whether to plot inferred CN or changepoint heatmap')

    acc_cmap = get_acc_cmap()
    plot_clustered_cell_cn_matrix(ax[1], df, 'rt_state_diff', cluster_field_name='clone_id', secondary_field_name='true_t', cmap=acc_cmap)
    ax[1].set_title('Replication accuracy: {}'.format(round(accuracy, 3)))
    

    fig.savefig(argv.output_rt_accuracy, bbox_inches='tight')


def plot_frac_rt_distributions(df, argv):
    df_frac = df[['cell_id', 'true_t', argv.frac_rt_col]].drop_duplicates().reset_index(drop=True)
    df_frac2 = df_frac.rename(columns={'true_t': 'true', argv.frac_rt_col: 'inferred'})
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
    set1 = df_frac['true_t'].values
    set2 = df_frac[argv.frac_rt_col].values
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

    # 1 is false positive, 0 is accurate, -1 is false negative
    df['rt_state_diff'] = df[argv.rep_col] - df['true_rep']

    # compute fraction of replicated bins per cells 
    if argv.frac_rt_col not in df.columns:
        df = compute_cell_frac(df, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.rep_col)

    # create separate plots
    plot_true_vs_inferred_rt_state(df, argv)
    plot_rt_accuracy(df, argv)
    plot_frac_rt_distributions(df, argv)


if __name__ == '__main__':
    main()
