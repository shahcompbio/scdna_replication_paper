from argparse import ArgumentParser
import pandas as pd
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


def main():
    argv = get_args()

    # load all datasets into one DataFrame
    df = load_data(argv)

    # create violinplot with t-tests for significance
    fig, ax = plt.subplots(1, 2, figsize=(8,4), tight_layout=True)
    ax = ax.flatten()

    # plot violins with p-values
    plot_clone_rt_diff_vs_cna_types(df, ax[0])
    ax[0].set_title('RT shifts at subclonal CNAs - all tumors')
    ax[0].set_xlabel('Subclonal CNA type')
    ax[0].set_ylabel('Clone RT relative to reference\n<--later | earlier-->')
    # draw a dashed horizontal line at y=0
    ax[0].axhline(y=0, linestyle='--', color='black')

    # plot histograms
    sns.histplot(
        data=df, x='clone_rt_diff', hue='clone_cna_type', ax=ax[1], 
        palette=get_cna_cmap(), alpha=0.5, hue_order=['loss', 'unaltered', 'gain'],
        common_norm=False, stat='density', bins=20
    )
    # draw a dashed vertical line at x=0
    ax[1].axvline(x=0, linestyle='--', color='black')
    ax[1].set_title('RT shifts at subclonal CNAs - all tumors')
    ax[1].set_xlabel('Clone RT relative to reference\n<--later | earlier-->')
    ax[1].set_ylabel('Density')

    fig.savefig(argv.plot, bbox_inches='tight', dpi=300)

    # save the merged table
    df.to_csv(argv.table, index=False)


if __name__=='__main__':
    main()
