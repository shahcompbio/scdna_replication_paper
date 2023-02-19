import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser()
    parser.add_argument('input', type=str, help='Input tsv file with true and inferred phase information for each cell')
    parser.add_argument('plot1', type=str, help='Confusion matrix plot for all cells')
    parser.add_argument('plot2', type=str, help='Plot of phase accuracies in parameter sweep')
    parser.add_argument('plot3', type=str, help='Plot of phase accuracies for all datasets')
    parser.add_argument('plot4', type=str, help='Jointplot of true vs inferred fraction of replicated bins')
    return parser.parse_args()


def plot_confusion_matrix(df, argv):
    ''' Given a table of true and inferred phases, plot a confusion matrix with counts of each cell in each phase '''
    # Plot confusion matrix
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), tight_layout=True)
    sns.heatmap(pd.crosstab(df['true_phase'], df['PERT_phase']), annot=True, fmt='d', ax=ax, cmap='Blues')
    ax.set_xlabel('PERT phase')
    ax.set_ylabel('True phase')
    ax.set_title('All simulated cells')
    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def violins_with_pvals(df, x, y, hue, ax, box_pairs, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    """ Create a violinplot with p-values annotated. """
    sns.violinplot(data=df, x=x, y=y, ax=ax)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_cna_rate_phase_acc(df, ax, n=1, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the replications state accuracy vs cna rate at a fixed number of clones (n). '''
    x = "cell_cna_rate"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('num_clones=={}'.format(n))
    box_pairs = [
        ((0.02, 'PERT comp.'), (0.00, 'PERT comp.')),
        ((0.02, 'PERT comp.'), (0.05, 'PERT comp.')),
        ((0.00, 'PERT comp.'), (0.05, 'PERT comp.'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test, text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('{} clone(s)'.format(n))
    return ax


def plot_clone_effect_phase_acc(df, ax, rate=0.02, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the replication accuracy vs number of clones at a fixed cell cna rate. '''
    x = "num_clones"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('cell_cna_rate=={}'.format(rate)).query('num_clones<4')
    box_pairs = [
        ((1, 'PERT comp.'), (3, 'PERT comp.'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('Cell CNA rate {}'.format(rate))
    return ax


def plot_param_sweep(df, argv):
    """ Figure showing phase accuracy parameter sweep across num_clones and cell_cna_rate """
    fig, ax = plt.subplots(1, 5, figsize=(20, 4), tight_layout=True)
    ax = ax.flatten()

    # specify the method name to use as a hue
    # statannot requires hues as input
    df['method'] = 'PERT comp.'

    # top row shows accuracy at predicting replication states
    # show effect of varying cna rate at fixed number of clones
    plot_cna_rate_phase_acc(df, ax[0], n=1)
    plot_cna_rate_phase_acc(df, ax[1], n=3)
    # show effect of varying number of clones at fixed cna rates
    plot_clone_effect_phase_acc(df, ax[2], rate=0.0)
    plot_clone_effect_phase_acc(df, ax[3], rate=0.02)
    plot_clone_effect_phase_acc(df, ax[4], rate=0.05)

    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)


def plot_all_datatags(df, argv):
    """ First figure is a mixture of barplots, scatterplots, and violinplots showing the accuracies of the different methods. """
    fig, ax = plt.subplots(2, 2, figsize=(8, 8), tight_layout=True)

    # merge together the two supblots in the top left corner
    gs = ax[0, 0].get_gridspec()
    for a in ax[0, :]:
        a.remove()
    axbig_top_row = fig.add_subplot(gs[0, 0:2])

    # showing rep accuracies on the top row
    # barplots and scatterplots of cn and rep accuracies for each method, across different simulation params
    sns.barplot(data=df, x='datatag', y='phase_acc', hue='num_clones', ax=axbig_top_row)

    # scatterplots which use A and lambda as the size params
    sns.scatterplot(data=df, x='cell_cna_rate', y='phase_acc', hue='num_clones', size='alpha', ax=ax[1, 0])
    sns.scatterplot(data=df, x='cell_cna_rate', y='phase_acc', hue='num_clones', size='lamb', ax=ax[1, 1])

    fig.savefig(argv.plot3, bbox_inches='tight', dpi=300)


def plot_jointplot(df, argv):
    ''' Plot a jointplot of the PERT and true fraction of replicated bins per cell. Use the true phase as the hue. '''
    fig, ax = plt.subplots(1, 1, figsize=(8, 8), tight_layout=True)
    sns.jointplot(data=df, x='cell_frac_rep', y='true_t', hue='true_phase', ax=ax)
    fig.savefig(argv.plot4, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    df['PERT_phase'] = df['PERT_phase'].astype(str)
    df['true_phase'] = df['true_phase'].astype(str)
    df['cell_id'] = df['cell_id'].astype(str)

    # change column name lambda to lambd to avoid conflict with python keyword
    df.rename(columns={'lambda': 'lamb'}, inplace=True)

    # Plot confusion matrix
    plot_confusion_matrix(df, argv)

    # Plot parameter sweep
    plot_param_sweep(df, argv)

    # Plot all datatags
    plot_all_datatags(df, argv)

    # plot a jointplot of the PERT and true fraction of replicated bins per cell
    plot_jointplot(df, argv)


if __name__ == '__main__':
    main()
