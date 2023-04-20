import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_phase_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='tables with phase calls and cell features where each table is a unique dataset')
    p.add_argument('-p1', '--plot1', type=str, help='confusion matrix comparing the phase calls of all models')
    p.add_argument('-p2', '--plot2', type=str, help='violin plots of the features across different phase calls')

    return p.parse_args()


def plot_confusion_matrix(cn, argv):
    """ Plot a confusion matrix comparing laks_phase to pert_phase for each cell. """
    # create a figure
    fig, ax = plt.subplots(figsize=(4, 4))

    # plot a confusion matrix of cn where the rows are laks phase and the columns are pert phase
    sns.heatmap(pd.crosstab(cn['pert_phase'], cn['laks_phase']), annot=True, fmt='d', cmap='Blues', ax=ax)

    # rename the x and y axis labels
    ax.set_xlabel('Laks phase')
    ax.set_ylabel('PERT phase')
    ax.set_title(f'Unsorted hTERT & OV2295\n # cells')

    # save the figure
    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def plot_violin_plots(cn, argv):
    """ Use violin plots to show the distribution of features across different phase calls. """
    # list of columns to plot on the y-axis
    y_cols = [
        'corrected_breakpoints', 'ploidy', 'corrected_madn', 
        'rpm_auto_norm', 'gc_intercept', 'gc_slope',
    ]

    # create a dictionary that maps each y-column name to a y-axis label
    y_label_dict = {
        'corrected_breakpoints': 'HMMcopy breakpoints',
        'ploidy': 'HMMcopy ploidy',
        'corrected_madn': 'RPM median absolute deviation\nbetween neighboring bins (madn)',
        'rpm_auto_norm': 'RPM autocorrelation',
        'gc_intercept': 'GC bias intercept',
        'gc_slope': 'GC bias slope'
    }

    # create a matplotlib figure with 3 rows and 3 columns
    fig, axes = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
    ax = axes.flatten()

    phase_cmap = get_phase_cmap()

    # loop through each column in y_cols, plotting the violin plot on the corresponding axis
    for y_col in y_cols:
        sns.violinplot(data=cn, x='laks_phase', hue='pert_phase', y=y_col, ax=ax[y_cols.index(y_col)], palette=phase_cmap, linewidth=1)
        ax[y_cols.index(y_col)].set_xlabel('Laks et al phase')
        ax[y_cols.index(y_col)].set_ylabel(y_label_dict[y_col])
        ax[y_cols.index(y_col)].legend(loc='upper right', title='PERT phase')
        ax[y_cols.index(y_col)].set_title(f'Unsorted hTERT & OV2295')

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
