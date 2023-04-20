import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_phase_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('pyro_s_input', type=str, help='pyro model output for PERT s-phase cells')
    p.add_argument('pyro_g_input', type=str, help='pyro model output for PERT g1/2-phase cells')
    p.add_argument('pyro_lq_input', type=str, help='pyro model output for PERT low quality cells')
    p.add_argument('frac_rt_col', type=str, help='column name for the fraction of replicated bins in a cell')
    p.add_argument('pyro_rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('dataset', type=str)
    p.add_argument('plot1', type=str, help='figure showing histogram of ccc features based on cell cycle phase')
    p.add_argument('plot2', type=str, help='figure showing scatterplot of ccc features based on cell cycle phase')
    p.add_argument('plot3', type=str, help='figure showing confusion matrix of true vs PERT inferred cell cycle phase')

    return p.parse_args()


def get_cell_metrics(cn_pyro_s, cn_pyro_g, cn_pyro_lq, frac_rt_col='cell_frac_rep'):
    # get cell-level metrics with the updated cell cycle labels
    cell_metric_cols = [
        'cell_id', 'true_phase', 'PERT_phase', 'laks_phase', 'breakpoints', 
        'total_mapped_reads_hmmcopy', 'madn', 'lrs', frac_rt_col, 
        'corrected_madn', 'corrected_breakpoints'
    ]
    cell_metrics_pyro_s = cn_pyro_s[cell_metric_cols].drop_duplicates()
    cell_metrics_pyro_g = cn_pyro_g[cell_metric_cols].drop_duplicates()
    cell_metrics_pyro_lq = cn_pyro_lq[cell_metric_cols].drop_duplicates()

    cell_metrics = pd.concat([cell_metrics_pyro_s, cell_metrics_pyro_g], ignore_index=True)
    cell_metrics = pd.concat([cell_metrics, cell_metrics_pyro_lq], ignore_index=True)

    return cell_metrics


def plot_ccc_distributions(cell_metrics, argv):
    phase_cmap = get_phase_cmap()

    # plot histogram of each feature split by cell cycle state
    fig, ax = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    sns.histplot(data=cell_metrics, x='corrected_madn', hue='PERT_phase', multiple='stack', ax=ax[0], palette=phase_cmap)
    sns.histplot(data=cell_metrics, x='lrs', hue='PERT_phase', multiple='stack', ax=ax[1], palette=phase_cmap)
    sns.histplot(data=cell_metrics, x='corrected_breakpoints', hue='PERT_phase', multiple='stack', ax=ax[2], palette=phase_cmap)

    # set the dataset to the title of every subplot
    for a in ax:
        a.set_title(argv.dataset)

    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)

    # plot pairwise scatterplots of all three features colored by cell cycle cstate
    fig, ax = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(data=cell_metrics, x='corrected_madn', y='lrs', hue='PERT_phase', alpha=0.3, ax=ax[0], palette=phase_cmap)
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='lrs', hue='PERT_phase', alpha=0.3, ax=ax[1], palette=phase_cmap)
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='corrected_madn', hue='PERT_phase', alpha=0.3, ax=ax[2], palette=phase_cmap)
    sns.scatterplot(data=cell_metrics, x='corrected_madn', y='lrs', hue=argv.frac_rt_col, alpha=0.3, ax=ax[3])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='lrs', hue=argv.frac_rt_col, alpha=0.3, ax=ax[4])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='corrected_madn', hue=argv.frac_rt_col, alpha=0.3, ax=ax[5])

    # set the dataset to the title of every subplot
    for a in ax:
        a.set_title(argv.dataset)

    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)

    # plot a confusion matrix comparing true_phase, PERT_phase, and laks_phase
    fig, ax = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)

    # true (x-axis) vs PERT (y-axis) in ax[0]
    sns.heatmap(pd.crosstab(cell_metrics.PERT_phase, cell_metrics.true_phase), annot=True, fmt='d', ax=ax[0], cmap='Blues')
    ax[0].set_ylabel('PERT phase')
    ax[0].set_xlabel('True phase')
    ax[0].set_title(argv.dataset)

    # true (x-axis) vs laks (y-axis) in ax[1]
    sns.heatmap(pd.crosstab(cell_metrics.laks_phase, cell_metrics.true_phase), annot=True, fmt='d', ax=ax[1], cmap='Blues')
    ax[1].set_ylabel('Laks phase')
    ax[1].set_xlabel('True phase')
    ax[1].set_title(argv.dataset)

    # laks (x-axis) vs PERT (y-axis) in ax[2]
    sns.heatmap(pd.crosstab(cell_metrics.PERT_phase, cell_metrics.laks_phase), annot=True, fmt='d', ax=ax[2], cmap='Blues')
    ax[2].set_ylabel('PERT phase')
    ax[2].set_xlabel('Laks phase')
    ax[2].set_title(argv.dataset)

    fig.savefig(argv.plot3, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn_pyro_s = pd.read_csv(argv.pyro_s_input, sep='\t')
    cn_pyro_g = pd.read_csv(argv.pyro_g_input, sep='\t')
    cn_pyro_lq = pd.read_csv(argv.pyro_lq_input, sep='\t')

    # use ecf status to come up with good vs bad S-phase labels in PERT_phase column
    cell_metrics = get_cell_metrics(cn_pyro_s, cn_pyro_g, cn_pyro_lq, frac_rt_col=argv.frac_rt_col)

    # rename the 'LowQual' entries in column 'PERT_phase' to 'LQ'
    cell_metrics.loc[cell_metrics.PERT_phase == 'LowQual', 'PERT_phase'] = 'LQ'

    # generate and save plots
    plot_ccc_distributions(cell_metrics, argv)
