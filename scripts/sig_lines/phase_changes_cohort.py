import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='tables with phase calls and cell features where each table is a unique dataset')
    p.add_argument('-p1', '--plot1', type=str, help='confusion matrix comparing the phase calls of all models')
    p.add_argument('-p2', '--plot2', type=str, help='violin plots of the features across different phase calls')

    return p.parse_args()


def get_phase_cmap():
    ''' Global color map for cell cycle phases '''
    cmap = {
        'S': '#BA0021',  # red
        'G1/2': '#003263',  # dark blue
        'G1': '#003263',  # dark blue
        'G2': '#6CACE4',  # light blue
        'LQ': '#C4CED4'  # silver
    }
    return cmap


def plot_confusion_matrix(cn, argv):
    """ Plot a confusion matrix comparing laks_phase to pert_phase for each cell. """
    # create a figure
    fig, ax = plt.subplots(figsize=(4, 4))

    # plot a confusion matrix of cn where the rows are laks phase and the columns are pert phase
    sns.heatmap(pd.crosstab(cn['laks_phase'], cn['pert_phase']), annot=True, fmt='d', cmap='Blues', ax=ax)

    # rename the x and y axis labels
    ax.set_xlabel('PERT phase')
    ax.set_ylabel('Laks et al phase')
    ax.set_title(f'Signatures cell lines\n # cells')

    # save the figure
    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def plot_violin_plots(cn, argv):
    """ Use violin plots to show the distribution of features across different phase calls. """
    # list of columns to plot on the y-axis
    y_cols = [
        'corrected_breakpoints', 'corrected_madn', 'rpm_auto_norm',
        'rep_auto_norm', 'cell_frac_rep', 'quality',
        'gc_intercept', 'gc_slope', 'ploidy'
    ]

    # create a dictionary that maps each y-column name to a y-axis label
    y_label_dict = {
        'corrected_breakpoints': 'CN breakpoints',
        'corrected_madn': 'RPM median absolute deviation\nbetween neighboring bins',
        'rpm_auto_norm': 'RPM autocorrelation',
        'rep_auto_norm': 'PERT rep state autocorrelation',
        'cell_frac_rep': 'PERT fraction of replicated bins',
        'quality': 'Laks et al quality score',
        'gc_intercept': 'GC bias intercept',
        'gc_slope': 'GC bias slope',
        'ploidy': 'HMMcopy ploidy'
    }

    # create a matplotlib figure with 3 rows and 3 columns
    fig, axes = plt.subplots(3, 3, figsize=(12, 12), tight_layout=True)
    ax = axes.flatten()

    phase_cmap = get_phase_cmap()

    # loop through each column in y_cols, plotting the violin plot on the corresponding axis
    for y_col in y_cols:
        sns.violinplot(data=cn, x='laks_phase', hue='pert_phase', y=y_col, ax=ax[y_cols.index(y_col)], palette=phase_cmap, linewidth=1)
        ax[y_cols.index(y_col)].set_xlabel('Laks et al phase')
        ax[y_cols.index(y_col)].set_ylabel(y_label_dict[y_col])
        ax[y_cols.index(y_col)].legend(loc='upper right', title='PERT phase')

    # save the figure
    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)



def main():
    argv = get_args()

    # read in the tables
    cn = pd.concat([pd.read_csv(f, sep='\t') for f in argv.input], ignore_index=True)

    # plot the confusion matrix
    plot_confusion_matrix(cn, argv)

    # plot the violin plots
    plot_violin_plots(cn, argv)


if __name__ == '__main__':
    main()
