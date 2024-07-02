from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_phase_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('summary_input', help='table for the accruacy of each permuted dataset')
    p.add_argument('metrics_input', help='table containing per-cell metrics for cells in all permuted datasets')
    p.add_argument('summary_plots', help='plot showing the model accuracy for catching mislabeled cells')
    p.add_argument('ccc_plots', help='plot showing the ccc features for mislabeled cells')

    return p.parse_args()


def make_plots(legend_df, metrics_df, argv):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    phase_cmap = get_phase_cmap()

    # barplot of the fraction of G1/2-cells accurately removed out of all those with swapped FACS labels
    sns.barplot(data=legend_df, x='rate', y='accuracy', ax=ax[1], color='#C4CED4')
    ax[1].set_ylabel('Fraction of mislabeled\ncells detected by model')
    ax[1].set_xlabel('Fraction of FACS G1/2-phase cells\nmislabeled as S-phase')
    ax[1].set_title('PERT phase accuracy')

    # distribution of cell_frac_rep values based on the true FACS sorting states
    # copy true_cell_cycle_state to a new column named 'FACS phase'
    metrics_df['FACS phase'] = metrics_df['true_cell_cycle_state']
    sns.histplot(data=metrics_df.query("cell_cycle_state=='S'"), x='cell_frac_rep', hue='FACS phase', bins=20, multiple='stack', ax=ax[0], palette=phase_cmap)
    ax[0].set_title('Permutation of FACS labels for PERT initialization')
    ax[0].set_xlabel('Inferred fraction of replicated bins')
    ax[0].set_ylabel('# cells')

    fig.savefig(argv.summary_plots, bbox_inches='tight', dpi=300)

    # create a new column entitled 'PERT phase' that is 'G' when extreme_cell_frac is True and 'S' when extreme_cell_frac is False
    metrics_df['PERT phase'] = metrics_df['extreme_cell_frac'].apply(lambda x: 'G1/2' if x else 'S')

    # rename the columns in metrics df to make the plots more readable
    metrics_df.rename(columns={
        'is_s_phase_prob': 'Laks S-phase probability',
        'corrected_breakpoints': 'CN breakpoints',
        'corrected_madn': 'RPM median absolute deviation',
        'quality': 'Laks quality score',
        }, inplace=True)

    # scatterplots of ccc features for the the mislabled cells, colored by whether the model thinks the cell is replicating (extreme_cell_frac)
    fig, ax = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(data=metrics_df.query("permuted==True"), x='Laks S-phase probability', y='Laks quality score', hue='PERT phase', alpha=0.5, ax=ax[0], palette=phase_cmap)
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='Laks S-phase probability', y='RPM median absolute deviation', hue='PERT phase', alpha=0.5, ax=ax[1], palette=phase_cmap)
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='Laks S-phase probability', y='CN breakpoints', hue='PERT phase', alpha=0.5, ax=ax[2], palette=phase_cmap)
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='RPM median absolute deviation', y='Laks quality score', hue='PERT phase', alpha=0.5, ax=ax[3], palette=phase_cmap)
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='CN breakpoints', y='Laks quality score', hue='PERT phase', alpha=0.5, ax=ax[4], palette=phase_cmap)
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='RPM median absolute deviation', y='CN breakpoints', hue='PERT phase', alpha=0.5, ax=ax[5], palette=phase_cmap)

    for i in range(6):
        ax[i].set_title('Mislabeled FACS G1/2 cells')
        ax[i].legend(title='PERT phase')

    fig.savefig(argv.ccc_plots, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    legend_df = pd.read_csv(argv.summary_input, sep='\t')
    metrics_df = pd.read_csv(argv.metrics_input, sep='\t')

    make_plots(legend_df, metrics_df, argv)


if __name__ == '__main__':
    main()
