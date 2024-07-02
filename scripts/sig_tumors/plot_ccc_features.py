import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='all cells with ccc features calculated')
    p.add_argument('plot1', type=str, help='figure showing histogram of ccc features based on cell cycle phase and extreme_cell_frac status')
    p.add_argument('plot2', type=str, help='figure showing scatterplot of ccc features based on cell cycle phase and extreme_cell_frac status')

    return p.parse_args()


def get_cell_metrics(cn):
    # get cell-level metrics with the updated cell cycle labels
    cell_metric_cols = [
        'cell_id', 'library_id', 'clone_id', 'breakpoints', 
        'total_mapped_reads_hmmcopy', 'madn', 'lrs', 
        'corrected_madn', 'corrected_breakpoints',
        'is_s_phase', 'is_s_phase_prob', 'in_tree'
    ]
    cell_metrics = cn[cell_metric_cols].drop_duplicates().reset_index(drop=True)

    return cell_metrics


def plot_ccc_distributions(cell_metrics, argv):
    # plot histogram of each feature split by cell cycle state
    fig, ax = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
    ax = ax.flatten()

    sns.histplot(data=cell_metrics, x='corrected_madn', hue='is_s_phase', multiple='stack', ax=ax[0])
    sns.histplot(data=cell_metrics, x='lrs', hue='is_s_phase', multiple='stack', ax=ax[1])
    sns.histplot(data=cell_metrics, x='corrected_breakpoints', hue='is_s_phase', multiple='stack', ax=ax[2])
    sns.histplot(data=cell_metrics, x='corrected_madn', hue='in_tree', multiple='stack', ax=ax[3])
    sns.histplot(data=cell_metrics, x='lrs', hue='in_tree', multiple='stack', ax=ax[4])
    sns.histplot(data=cell_metrics, x='corrected_breakpoints', hue='in_tree', multiple='stack', ax=ax[5])

    fig.savefig(argv.plot1, bbox_inches='tight')

    # plot pairwise scatterplots of all three features colored by cell cycle cstate
    fig, ax = plt.subplots(3, 3, figsize=(12, 12), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(data=cell_metrics, x='corrected_madn', y='lrs', hue='is_s_phase', alpha=0.3, ax=ax[0])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='lrs', hue='is_s_phase', alpha=0.3, ax=ax[1])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='corrected_madn', hue='is_s_phase', alpha=0.3, ax=ax[2])
    sns.scatterplot(data=cell_metrics, x='corrected_madn', y='lrs', hue='in_tree', alpha=0.3, ax=ax[3])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='lrs', hue='in_tree', alpha=0.3, ax=ax[4])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='corrected_madn', hue='in_tree', alpha=0.3, ax=ax[5])
    sns.scatterplot(data=cell_metrics, x='corrected_madn', y='lrs', hue='library_id', alpha=0.3, ax=ax[6])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='lrs', hue='library_id', alpha=0.3, ax=ax[7])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='corrected_madn', hue='library_id', alpha=0.3, ax=ax[8])

    fig.savefig(argv.plot2, bbox_inches='tight')


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input)
   
    # use ecf status to come up with good vs bad S-phase labels in cell_cycle_state column
    cell_metrics = get_cell_metrics(cn)

    # generate and save plots
    plot_ccc_distributions(cell_metrics, argv)

