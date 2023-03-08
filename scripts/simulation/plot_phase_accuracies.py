import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_methods_cmap, get_phase_cmap


def get_args():
    parser = ArgumentParser()
    parser.add_argument('input', type=str, help='Input tsv file with true and inferred phase information for each cell')
    parser.add_argument('plot1', type=str, help='Confusion matrix plot for all cells')
    parser.add_argument('plot2', type=str, help='Plot of phase accuracies in parameter sweep')
    parser.add_argument('plot3', type=str, help='Jointplot of true vs inferred fraction of replicated bins')
    return parser.parse_args()


def plot_confusion_matrix(df, argv):
    ''' Given a table of true and inferred phases, plot a confusion matrix with counts of each cell in each phase '''
    # subset to just the rows with method=='PERT'
    df = df.query('method=="PERT"').query('lamb==0.75')
    # Plot confusion matrix
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), tight_layout=True)
    sns.heatmap(pd.crosstab(df['predicted_phase'], df['true_phase']), annot=True, fmt='d', ax=ax, cmap='Blues')
    ax.set_xlabel('True phase')
    ax.set_ylabel('PERT phase')
    ax.set_title('All simulated cells (lambda=0.75)')
    # center the y-axis tick labels
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, va='center')
    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def violins_with_pvals(df, x, y, hue, ax, box_pairs, test='t-test_ind', text_format='star', loc='inside', verbose=0, show_hue=True):
    """ Create a violinplot with p-values annotated. """
    if show_hue:
        palette = get_methods_cmap()
        sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, palette=palette, linewidth=1)
    else:
        sns.violinplot(data=df, x=x, y=y, ax=ax, linewidth=1)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test,
                        text_format=text_format, loc=loc, verbose=verbose)
    # set legend in top right corner
    ax.legend(loc='upper right')


def plot_cna_rate_phase_acc(df, ax, n=1, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs cna rate at a fixed number of clones (n). '''
    x = "cell_cna_rate"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.75').query('num_clones=={}'.format(n)).query('beta0==1.2')
    box_pairs = [
        ((0.02, 'PERT'), (0.02, 'laks')),
        ((0.00, 'PERT'), (0.00, 'laks')),
        ((0.05, 'PERT'), (0.05, 'laks'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test, text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('Sweep across cell CNA rate\n# of clones={}'.format(n))
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('Cell CNA rate')
    return ax


def plot_clone_effect_phase_acc(df, ax, rate=0.02, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs number of clones at a fixed cell cna rate. '''
    x = "num_clones"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.75').query('cell_cna_rate=={}'.format(rate)).query('num_clones<4').query('beta0==1.2')
    box_pairs = [
        ((1, 'PERT'), (1, 'laks')),
        ((3, 'PERT'), (3, 'laks'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('Sweep across # of clones\nCell CNA rate={}'.format(rate))
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('Number of clones')
    return ax


def plot_alpha_effect(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs alpha where the hue is cell cna rate. '''
    x = "alpha"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('lamb==0.75').query('beta0==1.2').query('num_clones<4').query('cell_cna_rate==0.02')
    box_pairs = [
        ((10.0, 'PERT'), (10.0, 'laks')),
        ((5.0, 'PERT'), (5.0, 'laks')),
        ((15.0, 'PERT'), (15.0, 'laks'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('Sweep across alpha')
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('alpha')


def plot_lambda_effect(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs lambda where the hue is cell cna rate. '''
    x = "lamb"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0').query('num_clones==1').query('alpha==10.0').query('beta0==1.2')
    box_pairs = [
        ((0.5, 'laks'), (0.5, 'PERT')),
        ((0.6, 'laks'), (0.6, 'PERT')),
        ((0.75, 'laks'), (0.75, 'PERT')),
        ((0.9, 'laks'), (0.9, 'PERT')),
        ((0.99, 'laks'), (0.99, 'PERT'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_xlabel('lambda')
    ax.set_title('Sweep across lambda')
    ax.set_ylabel('Phase accuracy')


def plot_gc_bias_effect(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the phase accuracy vs GC bias coefficients in datasets with all other params fixed. '''
    x = "beta0"
    y = "phase_acc"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0.0').query('num_clones==1').query('alpha==10.0').query('lamb==0.75')
    box_pairs = [
        ((1.2, 'PERT'), (1.2, 'laks')),
        ((-1.2, 'PERT'), (-1.2, 'laks')),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    # sns.violinplot(x=x, y=y, data=temp_df, ax=ax)
    ax.set_title('Sweep across GC bias coefficients')
    ax.set_ylabel('Phase accuracy')
    ax.set_xlabel('Beta0 (GC bias slope)')
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
    sns.barplot(data=df, x='datatag', y='phase_acc', hue='method', ax=axbig_bottom_row, palette=get_methods_cmap(), saturation=1)
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
    ''' Plot a jointplot of the PERT and true fraction of replicated bins per cell. Use the phase class (TP, FP, TN, FN) as the hue. '''
    # subset to just the rows with method=='PERT'
    df = df.query('method=="PERT"').query('lamb==0.75')
    # phase_class_pal = {'TP': 'green', 'TN': 'blue', 'FP': 'red', 'FN': 'orange'}
    pal = get_phase_cmap()
    # create a JointGrid instance
    g = sns.JointGrid(data=df, x='true_cell_frac_rep', y='cell_frac_rep', hue='true_phase', palette=pal, height=4)
    # plot a scatterplot on the joint axes with alpha=0.2   
    # order the hues such that S is first, G1/2 is second
    g.plot_joint(sns.scatterplot, alpha=0.2, s=5)
    # plot a histogram of the x and y variables on the marginal axes and 20 bins
    g.plot_marginals(sns.histplot, kde=True, bins=20)
    # rename the axes
    g.set_axis_labels('True fraction of replicated bins', 'PERT inferred fraction of replicated bins')
    # rename the legend title to 'Phase class'
    g.ax_joint.legend(title='True phase')
    g.savefig(argv.plot3, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    df['PERT_phase'] = df['PERT_phase'].astype(str)
    df['true_phase'] = df['true_phase'].astype(str)
    df['laks_phase'] = df['laks_phase'].astype(str)
    df['cell_id'] = df['cell_id'].astype(str)

    # rename the lambda column to lambd to avoid conflict with python keyword
    df.rename(columns={'lambda': 'lamb'}, inplace=True)

    # melt `PERT_phase_acc` and `laks_phase_acc` into a single column `phase_acc`
    # and 'PERT_phase' and 'laks_phase' into a single column `pred_phase`
    # where the method is noted in a new column `method`
    df = pd.melt(df, 
        id_vars=[col for col in df.columns if not col.endswith('_phase_acc')],
        value_vars=['PERT_phase_acc', 'laks_phase_acc'], var_name='method', value_name='phase_acc')

    # rename the method column to remove the '_phase_acc' suffix
    df['method'] = df['method'].str.replace('_phase_acc', '')

    # create a new column named 'predicted_phase' that is the same as 'PERT_phase' if 'method' is 'PERT', 'laks_phase' if the method is 'laks'
    df['predicted_phase'] = np.where(df['method'] == 'PERT', df['PERT_phase'], df['laks_phase'])

    print(df[['cell_id', 'method', 'true_phase', 'predicted_phase']].head())

    # create a new column named 'phase_class' which says whether the 
    # given prediction is a true positive, false positive, true negative, or false negative
    df['phase_class'] = 'None'
    for i, row in df.iterrows():
        if row['predicted_phase'] == row['true_phase'] and row['true_phase'] == 'S':
            df.at[i, 'phase_class'] = 'TP'
        elif row['predicted_phase'] == row['true_phase'] and row['true_phase'] == 'G1/2':
            df.at[i, 'phase_class'] = 'TN'
        elif row['predicted_phase'] != row['true_phase'] and row['true_phase'] == 'S':
            df.at[i, 'phase_class'] = 'FN'
        elif row['predicted_phase'] != row['true_phase'] and row['true_phase'] == 'G1/2':
            df.at[i, 'phase_class'] = 'FP'

    print(df[['cell_id', 'method', 'true_phase', 'predicted_phase', 'phase_class']])
    print(df.phase_class.value_counts())

    # drop the pert_phase and laks_phase columns
    df.drop(columns=['PERT_phase', 'laks_phase'], inplace=True)
    
    # fill missing true_cell_frac_rep values with 0.0 as these are G1/2 cells with no replicated bins
    df['true_cell_frac_rep'].fillna(0.0, inplace=True)

    # rename 'LowQual' entries in the 'predicted_phase' column to 'LQ'
    df['predicted_phase'] = df['predicted_phase'].str.replace('LowQual', 'LQ')

    # Plot confusion matrix
    plot_confusion_matrix(df, argv)

    # Plot parameter sweep
    plot_param_sweep(df, argv)

    # plot a jointplot of the PERT and true fraction of replicated bins per cell
    plot_jointplot(df, argv)


if __name__ == '__main__':
    main()
