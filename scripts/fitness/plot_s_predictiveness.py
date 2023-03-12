import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LinearRegression
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_clone_cmap, get_rx_cmap



def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', nargs='+', help='list of cell cycle clone counts for each sample')
    p.add_argument('-p1', '--plot1', help='S-predictiveness plot for each sample')
    p.add_argument('-p2', '--plot2', help='S-predictiveness plot for all samples, split by Rx status')
    p.add_argument('-ot', '--output_tsv', help='table of all the data used to make the plots')

    return p.parse_args()


def timepoint_to_int(df):
    """Converts the timepoint column to an integer"""
    # get the timepoint column
    timepoints = df['timepoint'].values
    # convert the timepoints to integers
    timepoints = [int(x[1:]) for x in timepoints]
    # add the timepoints to the dataframe
    df['timepoint_int'] = timepoints
    # return the dataframe
    return df


def sort_timepoints(df):
    """Sort the dataframe according to clone_id and timepoint_int"""
    # if timepoint_int is not in the dataframe, create such a column
    df = timepoint_to_int(df)
    # sort the dataframe
    df = df.sort_values(by=['clone_id', 'timepoint_int'])
    # return the dataframe
    return df


def filter_rows(df, num_cells=10):
    """Filters out rows that do not have a value for instantaneous_s or have few cells"""
    df = df.loc[df['instantaneous_s'].notna()]
    df = df.loc[df['num_cells_g'] > num_cells]
    return df


def add_instantaneous_s_and_enrichment(df):
    """Adds a column to the dataframe that contains the observed clone shift in G1/2 population for each clone at each timepoint"""
    # compute the enrichment or depletion for S-phase cells at a given timepoint
    df['clone_frac_diff'] = df['clone_frac_s'] - df['clone_frac_g']

    times = sorted(df.timepoint_int.unique())
    clones = sorted(df.clone_id.unique())

    # this column is a proxy for that clone's instatneous selection coefficient
    df['instantaneous_s'] = np.nan

    # find difference in a clone's number/fraction of cells between two adjacent timepoints
    for t in range(len(times)-1):
        for c in clones:
            t0 = times[t]
            t1 = times[t+1]
            # find the row that corresponds to this clone & time
            row_t0 = df.loc[(df['clone_id']==c) & (df['timepoint_int']==t0)]
            row_t1 = df.loc[(df['clone_id']==c) & (df['timepoint_int']==t1)]
            # find the difference in G1/2-phase fractions between t0 and t1
            frac_diff = row_t1['clone_frac_g'].values[0] - row_t0['clone_frac_g'].values[0]

            # add frac_diff to the dataframe at the appropriate row
            df.loc[(df['clone_id']==c) & (df['timepoint_int']==t0), 'instantaneous_s'] = frac_diff
    return df


def plot_s_predictiveness(df, ax=None, title=None):
    """Plots the observed clone shift in G1/2 population vs. the clone's S-phase enrichment/depletion"""
    # if ax is None, create a new figure
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))

    # fit a regression line to the data
    sns.regplot(y='instantaneous_s', x='clone_frac_diff', data=df, ax=ax, scatter=False, color='black')

    # create a seaborn scatterplot comparing the observed clone shift in G1/2 population to the clone's S-phase enrichment/depletion
    sns.scatterplot(y='instantaneous_s', x='clone_frac_diff', data=df, hue='clone_id', style='timepoint', palette=get_clone_cmap(), ax=ax)
    # set the y-axis label
    ax.set_ylabel('Clone\n<-contraction | expansion->')
    # set the x-axis label
    ax.set_xlabel('S-phase\n<-depletion | enrichment->')
    # set the title
    if title is not None:
        ax.set_title(title)

    # expand the x axis limits to be slightly larger than the data
    ax.set_xlim(left=ax.get_xlim()[0] - 0.05, right=ax.get_xlim()[1] + 0.05)



# create a new plotting function that colors the data points by sample_id and uses different markers for the treatment status
def plot_s_predictiveness_cisplatin_combined(df, argv, title=None):
    """Plots the observed clone shift in G1/2 population vs. the clone's S-phase enrichment/depletion"""
    # fit a regression line to the data
    sns.lmplot(y='instantaneous_s', x='clone_frac_diff', data=df, scatter=True, hue='cisplatin', palette=get_rx_cmap())

    # set the y-axis label
    plt.ylabel('Clone\n<-contraction | expansion->')
    # set the x-axis label
    plt.xlabel('S-phase\n<-depletion | enrichment->')
    # set the title
    if title is not None:
        plt.title(title)

    plt.savefig(argv.plot2, dpi=300, bbox_inches='tight')


def main():
    argv = get_args()

    # create a list to store the dataframes for each sample
    df_pdx_combined = []

    # create a panel of figures to plot the S-predictiveness for each sample
    nrows = 8
    ncols = 1
    fig, ax = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows), tight_layout=True)
    ax = ax.flatten()

    # loop through the samples and plot their S-predictiveness in separate subpanels
    for i, path in enumerate(argv.input):
        # read in the data
        temp_df = pd.read_csv(path, sep='\t')
        # get the sample name, strip cisplatin suffix from SA535 samples
        sample = path.split('/')[2].replace('_CISPLATIN_Combined', '')

        # process the data and plot the s-predictiveness for this sample
        temp_df = sort_timepoints(temp_df)
        temp_df = add_instantaneous_s_and_enrichment(temp_df)
        temp_df = filter_rows(temp_df)
        plot_s_predictiveness(temp_df, ax=ax[i], title=sample)

        # save the sample_id in a new column
        temp_df['sample_id'] = sample

        # if 'U' is in the sample name, add a column that indicates that the sample is 'Rx-', otherwise it is 'Rx+'
        if 'U' in sample:
            temp_df['cisplatin'] = 'Rx-'
        else:
            temp_df['cisplatin'] = 'Rx+'
        
        # add the dataframe to the list
        df_pdx_combined.append(temp_df)
    
    # save the first figure
    fig.savefig(argv.plot1, dpi=300, bbox_inches='tight')

    # concatenate the dataframes into a single dataframe
    df_pdx_combined = pd.concat(df_pdx_combined, ignore_index=True)

    # plot the combined data and save the figure
    plot_s_predictiveness_cisplatin_combined(df_pdx_combined, argv, title='TNBC PDXs')

    # save the dataframe for posterity
    df_pdx_combined.to_csv(argv.output_tsv, sep='\t', index=False)



if __name__ == '__main__':
    main()
