import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='Table of the number of cells per cell cycle phase and clone')
    p.add_argument('output_clones', type=str, help='figure showing the evolution of clones over time for each phase')
    p.add_argument('output_total', type=str, help='figure showing total number of cells at each timepoint')

    return p.parse_args()


def plot_clone_stackplots(df, argv):
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True, sharey=True)

    clone_frac_g_vs_time = np.zeros((len(df.clone_id.unique()), len(df.timepoint.unique())))
    clone_frac_s_vs_time = np.zeros((len(df.clone_id.unique()), len(df.timepoint.unique())))
    clone_legend = []
    timepoints = df.timepoint.unique()
    i = 0
    for clone_id, chunk in df.groupby('clone_id'):
        chunk.sort_values(by='timepoint', inplace=True)
        clone_frac_g_vs_time[i] = chunk['clone_frac_g'].values
        clone_frac_s_vs_time[i] = chunk['clone_frac_s'].values
        clone_legend.append(clone_id)
        i += 1


    ax[0].stackplot(timepoints, clone_frac_g_vs_time, labels=clone_legend)
    ax[1].stackplot(timepoints, clone_frac_s_vs_time, labels=clone_legend)
    ax[0].set_ylabel('Clone fraction')
    ax[0].set_title('G1/2-phase')
    ax[1].set_title('S-phase')
    ax[1].legend(title='Clone ID')

    fig.savefig(argv.output_clones, bbox_inches='tight')


def plot_total_cell_count_stackplot(df, argv):
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    timepoints = df.timepoint.unique()
    total_counts = np.zeros((2, len(timepoints)))
    phase_legend = ['S', 'G1/2']

    i = 0
    for time, chunk in df.groupby('timepoint'):
        total_counts[0, i] = sum(chunk['num_cells_s'].values)
        total_counts[1, i] = sum(chunk['num_cells_g'].values)
        i += 1

    ax.stackplot(timepoints, total_counts, labels=phase_legend)
    ax.set_ylabel('# cells')
    ax.set_title('Total cell count per timepoint')
    ax.legend(title='Phase')

    fig.savefig(argv.output_total, bbox_inches='tight')


def main():
    argv = get_args()
    # load long-form dataframes from different cell cycle phases
    df = pd.read_csv(argv.input, sep='\t')

    # plot clone fractions for each phase & timepoint in the form of a stackplot
    plot_clone_stackplots(df, argv)

    # plot the total number of cells in each cell cycle phase for each timepoint
    # this should be helpful for identifying timepoints that should be removed
    plot_total_cell_count_stackplot(df, argv)


if __name__ == '__main__':
    main()
