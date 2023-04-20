import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_clone_cmap, get_htert_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('SA039_input', type=str, help='cell cycle clone counts from SA039')
    p.add_argument('SA906a_input', type=str, help='cell cycle clone counts from SA906a')
    p.add_argument('SA906b_input', type=str, help='cell cycle clone counts from SA906b')
    p.add_argument('output_tsv', type=str, help='Table of the number of cells per cell cycle phase and clone across all datasets after filtering')
    p.add_argument('output_png', type=str, help='figure comparing the instantaneous fitness to cell cycle enrichment for each clone across all datasets')

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


def add_instantaneous_s_and_enrichment(df):
    """Adds a column to the dataframe that contains the instantaneous selection coefficient for each clone at each timepoint"""
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


def filter_rows(df, num_cells=10):
    """Filters out rows that do not have a value for instantaneous_s or have few cells"""
    df = df.loc[df['instantaneous_s'].notna()]
    df = df.loc[df['num_cells_g'] > num_cells]
    return df


def plot_s_predictiveness(df, ax=None, title=None):
    """Plots the instantaneous selection coefficient vs. the clone's S-phase enrichment/depletion"""
    # if ax is None, create a new figure
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    # fit a regression line to the data
    sns.regplot(y='instantaneous_s', x='clone_frac_diff', data=df, ax=ax, scatter=False, color='black')

    # create a seaborn scatterplot comparing the instantaneous selection coefficient to the clone's S-phase enrichment/depletion
    sns.scatterplot(y='instantaneous_s', x='clone_frac_diff', data=df, hue='clone_id', style='timepoint', ax=ax, palette=get_clone_cmap())
    # set the y-axis label
    ax.set_ylabel('Clone\n<-contraction | expansion->')
    # set the x-axis label
    ax.set_xlabel('S-phase\n<-depletion | enrichment->')
    # set the title
    if title is not None:
        ax.set_title(title)

    # expand the x axis limits to be slightly larger than the data
    ax.set_xlim(left=ax.get_xlim()[0] - 0.05, right=ax.get_xlim()[1] + 0.05)


# create a new plotting function that colors the data points by sample_id
def plot_s_predictiveness_combined(df, ax=None, title=None):
    """Plots the instantaneous selection coefficient vs. the clone's S-phase enrichment/depletion"""
    # if ax is None, create a new figure
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))

    # fit a regression line to the data
    sns.regplot(y='instantaneous_s', x='clone_frac_diff', data=df, ax=ax, scatter=False, color='black')

    # create a seaborn scatterplot comparing the instantaneous selection coefficient to the clone's S-phase enrichment/depletion
    sns.scatterplot(y='instantaneous_s', x='clone_frac_diff', data=df, hue='sample_id', alpha=0.5, ax=ax, palette=get_htert_cmap())
    # set the y-axis label
    ax.set_ylabel('Clone\n<-contraction | expansion->')
    # set the x-axis label
    ax.set_xlabel('S-phase\n<-depletion | enrichment->')
    # set the title
    if title is not None:
        ax.set_title(title)

    # expand the x axis limits to be slightly larger than the data
    ax.set_xlim(left=ax.get_xlim()[0] - 0.05, right=ax.get_xlim()[1] + 0.05)


def main():
    argv = get_args()

    # read in the data
    df_SA039 = pd.read_csv(argv.SA039_input, sep='\t')
    df_SA906a = pd.read_csv(argv.SA906a_input, sep='\t')
    df_SA906b = pd.read_csv(argv.SA906b_input, sep='\t')

    # analyze the SA906a data
    # sort the timepoints chronologically
    df_SA906a = sort_timepoints(df_SA906a)
    # remove SA906a timepoints that are before X25 as they are prior to CRISPR perturbation and splitting into two lines
    df_SA906a = df_SA906a.loc[df_SA906a['timepoint_int'] >= 25]
    # add the instantaneous selection coefficient and S-phase enrichment/depletion columns
    df_SA906a = add_instantaneous_s_and_enrichment(df_SA906a)
    # remove rows that do not have a value for instantaneous_s or have few cells
    df_SA906a = filter_rows(df_SA906a)

    # analyze the SA906b data
    # sort the timepoints chronologically
    df_SA906b = sort_timepoints(df_SA906b)
    # remove SA906b timepoints that are before X30 as they are prior splitting into two lines
    df_SA906b = df_SA906b.loc[df_SA906b['timepoint_int'] >= 30]
    # add the instantaneous selection coefficient and S-phase enrichment/depletion columns
    df_SA906b = add_instantaneous_s_and_enrichment(df_SA906b)
    # remove rows that do not have a value for instantaneous_s or have few cells
    df_SA906b = filter_rows(df_SA906b)

    # analyze the SA039 data
    # sort the timepoints chronologically
    df_SA039 = sort_timepoints(df_SA039)
    # add the instantaneous selection coefficient and S-phase enrichment/depletion columns
    df_SA039 = add_instantaneous_s_and_enrichment(df_SA039)
    # remove rows that do not have a value for instantaneous_s or have few cells
    df_SA039 = filter_rows(df_SA039)

    # create a new column in each dataframe that stores the sample_id
    df_SA039['sample_id'] = 'SA039'
    df_SA906a['sample_id'] = 'SA906a'
    df_SA906b['sample_id'] = 'SA906b'
    # combine the dataframes from SA039, SA906a, and SA906b
    df_combined = pd.concat([df_SA039, df_SA906a, df_SA906b])

    fig, ax = plt.subplots(2, 2, figsize=(10, 10), tight_layout=True)

    # plot the combined data in the top left subplot
    plot_s_predictiveness_combined(df_combined, ax=ax[0, 0], title='All hTERT time-series')
    # plot the SA039 data in the top right subplot
    plot_s_predictiveness(df_SA039, ax=ax[0, 1], title='SA039 (WT)')
    # plot the SA906a data in the bottom left subplot
    plot_s_predictiveness(df_SA906a, ax=ax[1, 0], title='SA906a (p53-/-)')
    # plot the SA906b data in the bottom right subplot
    plot_s_predictiveness(df_SA906b, ax=ax[1, 1], title='SA906b (p53-/-)')

    # save the figure
    fig.savefig(argv.output_png, dpi=300, bbox_inches='tight')

    # save df_combined as a tsv
    df_combined.to_csv(argv.output_tsv, sep='\t', index=False)



if __name__ == '__main__':
    main()
