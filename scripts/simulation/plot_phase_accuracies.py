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
    parser.add_argument('plot3', type=str, help='Jointplot of true vs inferred fraction of replicated bins')
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


def violins_with_pvals(df, x, y, hue, ax, box_pairs, test='t-test_ind', text_format='star', loc='inside', verbose=0, show_hue=True):
    """ Create a violinplot with p-values annotated. """
    if show_hue:
        sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax)
    else:
        sns.violinplot(data=df, x=x, y=y, ax=ax)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_cna_rate_phase_acc(df, ax, n=1, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs cna rate at a fixed number of clones (n). '''
    x = "cell_cna_rate"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('num_clones=={}'.format(n))
    box_pairs = [
        ((0.02, 'PERT comp.'), (0.00, 'PERT comp.')),
        ((0.02, 'PERT comp.'), (0.05, 'PERT comp.')),
        ((0.00, 'PERT comp.'), (0.05, 'PERT comp.'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test, text_format=text_format, loc=loc, verbose=verbose, show_hue=False)
    ax.set_title('{} clone(s)'.format(n))
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('Cell CNA rate')
    return ax


def plot_clone_effect_phase_acc(df, ax, rate=0.02, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs number of clones at a fixed cell cna rate. '''
    x = "num_clones"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('cell_cna_rate=={}'.format(rate)).query('num_clones<4')
    box_pairs = [
        ((1, 'PERT comp.'), (3, 'PERT comp.'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, show_hue=False)
    ax.set_title('Cell CNA rate {}'.format(rate))
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('Number of clones')
    return ax


def plot_alpha_effect(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs alpha where the hue is cell cna rate. '''
    x = "cell_cna_rate"
    y = "phase_acc"
    hue = "alpha"
    temp_df = df.copy()
    box_pairs = [
        ((0.0, 5), (0.0, 10)),
        ((0.0, 5), (0.0, 15)),
        ((0.0, 10), (0.0, 15)),
        ((0.02, 5), (0.02, 10)),
        ((0.02, 5), (0.02, 15)),
        ((0.02, 10), (0.02, 15)),
        ((0.05, 5), (0.05, 10)),
        ((0.05, 5), (0.05, 15)),
        ((0.05, 10), (0.05, 15))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('Diploid & Polyclonal')
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('Cell CNA rate')


def plot_lambda_effect(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs lambda where the hue is cell cna rate. '''
    x = "lamb"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0').query('num_clones==1').query('alpha==10.0')
    box_pairs = [
        ((0.1, 'PERT comp.'), (0.5, 'PERT comp.')),
        ((0.1, 'PERT comp.'), (0.7, 'PERT comp.')),
        ((0.1, 'PERT comp.'), (0.85, 'PERT comp.')),
        ((0.1, 'PERT comp.'), (0.99, 'PERT comp.')),
        ((0.5, 'PERT comp.'), (0.7, 'PERT comp.')),
        ((0.5, 'PERT comp.'), (0.85, 'PERT comp.')),
        ((0.5, 'PERT comp.'), (0.99, 'PERT comp.')),
        ((0.7, 'PERT comp.'), (0.85, 'PERT comp.')),
        ((0.7, 'PERT comp.'), (0.99, 'PERT comp.')),
        ((0.85, 'PERT comp.'), (0.99, 'PERT comp.'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, show_hue=False)
    ax.set_xlabel('lambda')
    ax.set_title('Diploid, cell CNA rate 0, alpha 10')
    ax.set_ylabel('Phase accuracy')


def plot_gc_bias_effect(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs GC bias coefficients in datasets with all other params fixed. '''
    # note that the GC bias coefficients are 1.2,0.0 for all datasets except for D4 which is 0.1,0.2,-1,1,-0.25
    df['gc_bias'] = '1st order'
    df.loc[df['datatag'] == 'D4', 'gc_bias'] = '4th order'
    x = "gc_bias"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0.0').query('num_clones==1').query('alpha==10.0').query('lamb==0.7')
    box_pairs = [
        (('1st order', 'PERT comp.'), ('4th order', 'PERT comp.'))
    ]
    # violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
    #                    text_format=text_format, loc=loc, verbose=verbose, show_hue=False)
    sns.violinplot(x=x, y=y, data=temp_df, ax=ax)
    ax.set_title('Diploid, cell CNA rate 0.0\nalpha 10, lambda 0.7')
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('GC bias polynomial coefficients')
    # ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')


def plot_param_sweep(df, argv):
    """ Figure showing phase accuracy parameter sweep across num_clones and cell_cna_rate """
    fig, ax = plt.subplots(2, 5, figsize=(20, 8), tight_layout=True)

    # top row shows accuracy at predicting replication states
    # show effect of varying cna rate at fixed number of clones
    plot_cna_rate_phase_acc(df, ax[0, 0], n=1)
    plot_cna_rate_phase_acc(df, ax[0, 1], n=3)
    # show effect of varying number of clones at fixed cna rates
    plot_clone_effect_phase_acc(df, ax[0, 2], rate=0.0)
    plot_clone_effect_phase_acc(df, ax[0, 3], rate=0.02)
    plot_clone_effect_phase_acc(df, ax[0, 4], rate=0.05)

    # merge together the two supblots in the bottom left corner
    gs = ax[1, 0].get_gridspec()
    for a in ax[1, :2]:
        a.remove()
    axbig_bottom_row = fig.add_subplot(gs[1, 0:2])

    # barplots of phase accuracies for all simulated datasets
    sns.barplot(data=df, x='datatag', y='phase_acc', hue='num_clones', ax=axbig_bottom_row)
    axbig_bottom_row.set_ylabel('Phase accuracy')
    axbig_bottom_row.set_title('All simulated datasets')

    # plot the effect of varying alpha
    plot_alpha_effect(df, ax[1, 2])

    # plot the effect of varying lambda
    plot_lambda_effect(df, ax[1, 3])
    
    # plot the effect of varying GC bias
    plot_gc_bias_effect(df, ax[1, 4])

    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)



def plot_jointplot(df, argv):
    ''' Plot a jointplot of the PERT and true fraction of replicated bins per cell. Use the true phase as the hue. '''
    # create a JointGrid instance
    g = sns.JointGrid(data=df, x='cell_frac_rep', y='true_t', hue='true_phase', hue_order=['S', 'G1/2'])
    # plot a scatterplot on the joint axes with alpha=0.2   
    # order the hues such that S is first, G1/2 is second
    g.plot_joint(sns.scatterplot, alpha=0.2, s=5)
    # plot a histogram of the x and y variables on the marginal axes and 20 bins
    g.plot_marginals(sns.histplot, kde=True, bins=20)
    # rename the axes
    g.set_axis_labels('Inferred fraction of replicated bins', 'True fraction of replicated bins')
    g.savefig(argv.plot3, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    df['PERT_phase'] = df['PERT_phase'].astype(str)
    df['true_phase'] = df['true_phase'].astype(str)
    df['cell_id'] = df['cell_id'].astype(str)

    # specify the method name to use as a hue
    # statannot requires hues as input
    df['method'] = 'PERT comp.'

    # fill missing true_t values with 0.0 as these are G1/2 cells with no replicated bins
    df['true_t'].fillna(0.0, inplace=True)

    # change column name lambda to lambd to avoid conflict with python keyword
    df.rename(columns={'lambda': 'lamb'}, inplace=True)

    # Plot confusion matrix
    plot_confusion_matrix(df, argv)

    # Plot parameter sweep
    plot_param_sweep(df, argv)

    # plot a jointplot of the PERT and true fraction of replicated bins per cell
    plot_jointplot(df, argv)


if __name__ == '__main__':
    main()
