import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_rx_cmap


def get_args():
    parser = ArgumentParser()

    parser.add_argument('input', type=str, help='table of cell cycle counts for each clone-rx-dataset combination')
    parser.add_argument('output', type=str, help='output file name')

    return parser.parse_args()


def swarmplot_with_pvals(df, x, y, hue, ax, box_pairs, order=None, test='t-test_ind', text_format='star', loc='inside', verbose=0, palette=None):
    sns.swarmplot(data=df, x=x, y=y, hue=hue, ax=ax, order=order, palette=palette)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test, order=order,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_rx_spf(df, ax, x='rx_status', y='frac_cells_s', test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' 
    Plot the difference in pseudobulk clone S-phase fractions
    where the data is split by treatment status
    '''
    rx_cmap = get_rx_cmap()
    x = x
    y = y
    hue = None
    box_pairs = [
        ('untreated', 'treated'),
    ]
    order = ['untreated', 'treated']
    swarmplot_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose, palette=rx_cmap)
    return ax


def plot_spf_distributions(df, argv):

    fig, ax = plt.subplots(1, 1, figsize=(4, 4), tight_layout=True)

    # swarmplot with p-values
    plot_rx_spf(df, ax, x='rx_status', y='frac_cells_s')

    ax.set_ylabel('Clone S-phase fraction')
    ax.set_xlabel('')
    ax.set_title('Cell cycle distribution')

    fig.savefig(argv.output, bbox_inches='tight', dpi=300)


def filter_rows(df, threshold=10):
    ''' Filter the dataframe to only include rows if the clone & dataset has at least 10 cells in each rx_status'''
    # loop through each clone & dataset combination
    for (clone, dataset), clone_df in df.groupby(['clone_id', 'dataset']):
        # check to see taht there are both 'treated' and 'untreated' rx_status rows in clone_df
        if not {'treated', 'untreated'}.issubset(clone_df['rx_status'].unique()):
            # if not, drop the rows for this clone & dataset
            df = df.drop(clone_df.index)
            continue
        # find the number of S-phase cells in the treated and untreated groups
        treated_num_cells_s = clone_df.loc[clone_df['rx_status'] == 'treated', 'num_cells_s'].values[0]
        untreated_num_cells_s = clone_df.loc[clone_df['rx_status'] == 'untreated', 'num_cells_s'].values[0]
        # check if there are at least 10 cells in each rx_status
        if treated_num_cells_s < 10 or untreated_num_cells_s < 10:
            # if not, drop the rows for this clone & dataset
            df = df.drop(clone_df.index)
    
    return df


def simple_filter_rows(df, threshold=20):
    # filter all rows in df to have at least 20 cells between both phases
    df = df.loc[(df['num_cells_s'] + df['num_cells_g'] >= threshold)]
    return df


def main():
    argv = get_args()

    # load the data
    df = pd.read_csv(argv.input, sep='\t')

    # filter the data
    # df = filter_rows(df)
    df = simple_filter_rows(df)

    # plot
    plot_spf_distributions(df, argv)


if __name__ == '__main__':
    main()

