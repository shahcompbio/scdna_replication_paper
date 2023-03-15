from argparse import ArgumentParser
import numpy as np
import pandas as pd
from statannot import add_stat_annotation
import matplotlib.pyplot as plt
import seaborn as sns
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_rx_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='long-form scRT results every dataset')
    p.add_argument('--plot1', help='histogram of cell time in S-phase across the whole cohort')
    p.add_argument('--plot2', help='swarplot of cell time in S-phase, split by cisplatin treatment status and dataset')

    return p.parse_args()


def load_data(argv):
    """ Read in the scRT for each data cohort, only keeping the relevant per-cell columns. """
    df = []

    # loop through the input files
    for path in argv.input:
        # read in the tsv
        temp_df = pd.read_csv(path, sep='\t')
        # use the path name to extract the dataset name
        d = path.split('/')[2]
        # add the dataset name as a column
        temp_df['dataset'] = d
        # keep only the relevant columns
        temp_df = temp_df[['dataset', 'label', 'datasetname', 'cell_id', 'library_id', 'clone_id', 'cell_frac_rep', 'model_tau']].drop_duplicates().reset_index(drop=True)
        # add a 'cisplatin' column that is False if a 'U' appears in the datasetname of the cell, True otherwise
        temp_df['cisplatin'] = temp_df['datasetname'].apply(lambda x: 'untreated' if 'U' in x else 'treated')
        
        # append to the list
        df.append(temp_df)
    
    # concatenate the list of dataframes into one
    df = pd.concat(df, ignore_index=True)

    return df

def violinplot_with_pvals(df, x, y, hue, ax, box_pairs, order=None, test='t-test_ind', text_format='star', loc='inside', verbose=0, palette=None):
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, order=order, palette=palette)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test, order=order,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_rx_s_times(df, ax, x='dataset', y='cell_frac_rep', hue='cisplatin', test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' 
    Plot the difference in pseudobulk clone S-phase fractions
    where the data is split by treatment status
    '''
    rx_cmap = get_rx_cmap()
    x = x
    y = y
    hue = hue
    box_pairs = [
        (('SA535', 'untreated'), ('SA535', 'treated')),
        (('SA1035', 'untreated'), ('SA1035', 'treated')),
        (('SA609', 'untreated'), ('SA609', 'treated'))
    ]
    violinplot_with_pvals(df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, palette=rx_cmap)
    return ax


def main():
    argv = get_args()

    # load the data across the entire cohort
    df = load_data(argv)

    dataset_list = df['dataset'].unique()

    nrows = 2
    ncols = 3
    fig, ax = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows), tight_layout=True)

    # plot histograms of the cell S-phase times
    # top row is with common_norm=True, bottom row is with common_norm=False
    # plot a curve for every experiment (i.e. label)
    sns.histplot(data=df, x='cell_frac_rep', hue='label', ax=ax[0, 0])
    sns.histplot(data=df, x='cell_frac_rep', hue='label', common_norm=False, ax=ax[1, 0])

    # plot curves for on- vs off-cisplatin
    sns.histplot(data=df, x='cell_frac_rep', hue='cisplatin', ax=ax[0, 1])
    sns.histplot(data=df, x='cell_frac_rep', hue='cisplatin', common_norm=False, ax=ax[1, 1])

    # plot curves for each dataset
    sns.histplot(data=df, x='cell_frac_rep', hue='dataset', ax=ax[0, 2])
    sns.histplot(data=df, x='cell_frac_rep', hue='dataset', common_norm=False, ax=ax[1, 2])

    # set the x-axis labels and titles
    for i in range(nrows):
        for j in range(ncols):
            ax[i, j].set_xlabel('Inferred fraction of replicated loci')
            ax[i, j].set_title('Cell S-phase times')
    
    fig.savefig(argv.plot1, dpi=300, bbox_inches='tight')

    # plot the difference between on- vs off-cisplatin within each dataset
    fig, ax = plt.subplots(1, 1, figsize=(4, 4), tight_layout=True)
    plot_rx_s_times(df, ax)
    ax.set_xlabel('')
    ax.set_ylabel('Inferred fraction of replicated loci\n<-early | late->')
    ax.set_title('Cell S-phase times')
    ax.legend(loc='lower right')
    fig.savefig(argv.plot2, dpi=300, bbox_inches='tight')



if __name__ == '__main__':
    main()
