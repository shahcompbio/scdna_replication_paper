from argparse import ArgumentParser
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_cna_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='subclonal rt diff tables for every dataset')
    p.add_argument('--table', help='table of all the concatenated input tables')
    p.add_argument('--plot', help='Summary plots of subclonal RT diffs across all datasets')

    return p.parse_args()


def load_data(argv):
    # load subclonal rt diff table for each dataset
    df = []
    for path in argv.input:
        temp_df = pd.read_csv(path)
        df.append(temp_df)

    # concatenate into one df
    df = pd.concat(df, ignore_index=True)
    return df


def violins_with_pvals(df, x, y, hue, ax, box_pairs, order=None, test='t-test_ind', text_format='star', loc='inside', verbose=0, palette=None):
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, order=order, palette=palette)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test, order=order,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_clone_rt_diff_vs_cna_types(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the distribution of clone RT differences against the CNA type of that particular locus. '''
    cna_cmap = get_cna_cmap()
    x = "clone_cna_type"
    y = "clone_rt_diff"
    hue = None
    box_pairs = [
        ('loss', 'gain'),
        ('loss', 'unaltered'),
        ('unaltered', 'gain'),
    ]
    order = ['loss', 'unaltered', 'gain']
    violins_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose, palette=cna_cmap)
    return ax


def subclonal_rt_shift_histogram(df, ax, title_prefix):
    # show the histogram of clone_rt_diff with a hue of clone_cna_type
    sns.histplot(
        data=df, x='clone_rt_diff', hue='clone_cna_type', hue_order=['loss', 'unaltered', 'gain'], 
        palette=get_cna_cmap(), common_norm=False, stat='density', ax=ax,
        cbar_kws={'edgecolor': 'k', 'linewidth': 0}, kde=True
    )

    # compute the p-values between loss-unaltered and gain-unaltered
    # get the y-values for the loss-unaltered and gain-unaltered comparisons
    loss = df[df.clone_cna_type == 'loss'].clone_rt_diff
    unaltered = df[df.clone_cna_type == 'unaltered'].clone_rt_diff
    gain = df[df.clone_cna_type == 'gain'].clone_rt_diff
    # compute the p-value between loss-unaltered and gain-unaltered
    # multiply by 3 because we are doing 3 tests (bonferroni correction)
    loss_unaltered_pval = stats.ttest_ind(loss, unaltered)[1] * 3
    gain_unaltered_pval = stats.ttest_ind(gain, unaltered)[1] * 3
    gain_loss_pval = stats.ttest_ind(gain, loss)[1] * 3

    # annotate the p-values in the upper-left corner
    ax.text(0.05, 0.9, 'gain-unalt p={:.2e}'.format(gain_unaltered_pval), transform=ax.transAxes, ha='left', va='center')
    ax.text(0.05, 0.8, 'gain-loss p={:.2e}'.format(gain_loss_pval), transform=ax.transAxes, ha='left', va='center')
    ax.text(0.05, 0.7, 'loss-unalt p={:.2e}'.format(loss_unaltered_pval), transform=ax.transAxes, ha='left', va='center')

    # annotate the median of each distribution with a dashed vertical line
    # annotate the line with the median value as text with 2 decimal places
    loss_median = df[df.clone_cna_type=='loss'].clone_rt_diff.median()
    ax.axvline(x=loss_median, color=get_cna_cmap()['loss'], linestyle='--')
    # ax.text(x=loss_median, y=10, s='median={:.2f}'.format(loss_median), color=get_cna_cmap()['loss'], ha='left', va='center')
    unaltered_median = df[df.clone_cna_type=='unaltered'].clone_rt_diff.median()
    ax.axvline(x=unaltered_median, color=get_cna_cmap()['unaltered'], linestyle='--')
    # ax.text(x=unaltered_median, y=12, s='median={:.2f}'.format(unaltered_median), color=get_cna_cmap()['unaltered'], ha='left', va='center')
    gain_median = df[df.clone_cna_type=='gain'].clone_rt_diff.median()
    ax.axvline(x=gain_median, color=get_cna_cmap()['gain'], linestyle='--')
    # ax.text(x=gain_median, y=10, s='median={:.2f}'.format(gain_median), color=get_cna_cmap()['gain'], ha='right', va='center')

    ax.set_title('{}: RT shifts at subclonal CNAs'.format(title_prefix))
    ax.set_xlabel('Clone RT relative to reference\n<--later | earlier-->')
    ax.set_xlim(-0.3, 0.3)
    # add xticks and xticklabels from -0.3 to 0.3 spaced by 0.1
    # round each value to 1 decimal place
    ax.set_xticks(np.arange(-0.3, 0.4, 0.1))
    ax.set_xticklabels(np.arange(-0.3, 0.4, 0.1).round(1))
    # remove the legend from the plot
    ax.legend().remove()


def main():
    argv = get_args()

    # load all datasets into one DataFrame
    df = load_data(argv)

    # create violinplot with t-tests for significance
    fig, ax = plt.subplots(1, 1, figsize=(6,4))
    # ax = plot_clone_rt_diff_vs_cna_types(df, ax)
    # ax.set_title('RT shifts at subclonal CNAs - all hTERTs')
    # ax.set_xlabel('Subclonal CNA type')
    # ax.set_ylabel('Clone RT relative to reference\n<--later | earlier-->')
    subclonal_rt_shift_histogram(df, ax, 'Polyclonal simulated data')
    fig.savefig(argv.plot, bbox_inches='tight', dpi=300)

    # save the merged table
    df.to_csv(argv.table, sep='\t', index=False)


if __name__=='__main__':
    main()
