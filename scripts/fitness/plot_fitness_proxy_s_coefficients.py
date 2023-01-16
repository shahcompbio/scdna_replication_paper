import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('SA1035_treated_input', type=str, help='Table of the number of cells per cell cycle phase and clone for the SA1035 treated sample')
    p.add_argument('SA1035_untreated_input', type=str, help='Table of the number of cells per cell cycle phase and clone for the SA1035 untreated sample')
    p.add_argument('SA535_treated_input', type=str, help='Table of the number of cells per cell cycle phase and clone for the SA535 treated sample')
    p.add_argument('SA535_untreated_input', type=str, help='Table of the number of cells per cell cycle phase and clone for the SA535 untreated sample')
    p.add_argument('output', type=str, help='figure showing the relative fitness coefficients for each clone between the treated and untreated samples')

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


def fill_in_missing_clones(df, clone_list=None):
    """ If there are clones present at some timepoints but not others, fill in the missing clones with 0s """
    # find the clone list if not specified
    if clone_list is None:
        clone_list = df.clone_id.unique()
    # loop through the clones and timepoints
    for clone in clone_list:
        for timepoint in df.timepoint.unique():
            # if the clone is not present at the timepoint, add a row with 0s
            if clone not in df[df.timepoint == timepoint].clone_id.values:
                new_line = pd.DataFrame({
                    'clone_id': [clone], 'timepoint': [timepoint], 'num_cells_s': [0], 'num_cells_g': [0],
                    'clone_frac_s': [0], 'clone_frac_g': [0], 'positive_p': [0], 'negative_p': [0],
                    'positive_p_adj': [0], 'negative_p_adj': [0], 'timepoint_int': [int(timepoint[1:])]
                })
                # concatenate the new line to the dataframe
                df = pd.concat([df, new_line], ignore_index=True)
    return df


def calc_proxy_s(df):
    """Calculate the proxy s coefficient for each clone use the G1/2-phase fractions of the first and last timepoints"""
    # subset the data to only include the first and last timepoints
    df = sort_timepoints(df)
    times = sorted(df.timepoint_int.unique())
    df2 = df[df.timepoint_int.isin([times[0], times[-1]])]
    # print(df2)
    
    # subtract clone_frac_g for each clone between the first and last timepoints
    proxy_s_df = []
    for clone_id, chunk in df2.groupby('clone_id'):
        print(chunk)
        temp_diff = chunk.clone_frac_g.values[1] - chunk.clone_frac_g.values[0]
        proxy_s_df.append(pd.DataFrame({'clone_id': [clone_id], 'proxy_s': [temp_diff]}))
    proxy_s_df = pd.concat(proxy_s_df, ignore_index=True)
    
    return proxy_s_df


def plot_proxy_s_coefficients(proxy_s_df_SA1035, proxy_s_df_SA535, argv):
    """
    Given the clone proxy s coefficients for sample, plot show the split between on and off treatment using a pointplot.
    """
    # create a 2x2 grid of plots
    fig, ax = plt.subplots(2, 2, figsize=(10, 10), tight_layout=True)
    ax = ax.flatten()

    # plot the proxy_s_rank in the top row
    sns.pointplot(data=proxy_s_df_SA1035, x='treatment', y='proxy_s_rank', hue='clone_id', ax=ax[0])
    ax[0].set_ylabel('proxy s rank\n<- low fitness | high fitness ->')
    ax[0].set_title('SA1035')
    sns.pointplot(data=proxy_s_df_SA535, x='treatment', y='proxy_s_rank', hue='clone_id', ax=ax[1])
    ax[1].set_ylabel('proxy s rank\n<- low fitness | high fitness ->')
    ax[1].set_title('SA535')

    # plot the raw proxy_s in the bottom row
    sns.pointplot(data=proxy_s_df_SA1035, x='treatment', y='proxy_s', hue='clone_id', ax=ax[2])
    ax[2].set_ylabel('proxy s\n<- low fitness | high fitness ->')
    ax[2].set_title('SA1035')
    sns.pointplot(data=proxy_s_df_SA535, x='treatment', y='proxy_s', hue='clone_id', ax=ax[3])
    ax[3].set_ylabel('proxy s\n<- low fitness | high fitness ->')
    ax[3].set_title('SA535')

    # save the figure
    plt.savefig(argv.output, dpi=300, bbox_inches='tight')


def main():
    argv = get_args()

    # load the cell cycle clone counts for all the samples
    df_SA1035U = pd.read_csv(argv.SA1035_untreated_input, sep='\t')
    df_SA1035T = pd.read_csv(argv.SA1035_treated_input, sep='\t')
    df_SA535U = pd.read_csv(argv.SA535_untreated_input, sep='\t')
    df_SA535T = pd.read_csv(argv.SA535_treated_input, sep='\t')

    # sort timepoints based on the timepoint_int column
    df_SA1035U = sort_timepoints(df_SA1035U)
    df_SA1035T = sort_timepoints(df_SA1035T)
    df_SA535U = sort_timepoints(df_SA535U)
    df_SA535T = sort_timepoints(df_SA535T)

    # add the pretretament timepoint for SA1035
    # find the earliest timepoint in the untreated sample
    earliest_timepoint_SA1035 = sorted(df_SA1035U.timepoint_int.unique())[0]
    # add data from the earliest untreated timepiont to the treated dataframe
    df_SA1035T = pd.concat([df_SA1035U[df_SA1035U.timepoint_int == earliest_timepoint_SA1035], df_SA1035T], ignore_index=True)

    # add the pretreatement timepoint for SA535
    # find the earliest timepoint in the untreated sample
    earliest_timepoint_SA535 = sorted(df_SA535U.timepoint_int.unique())[0]
    # add data from the earliest untreated timepiont to the treated dataframe
    df_SA535T = pd.concat([df_SA535U[df_SA535U.timepoint_int == earliest_timepoint_SA535], df_SA535T], ignore_index=True)

    # fill in the missing clones for SA1035
    # find the set of clone_ids that appear in the union of the treated and untreated samples
    clone_list_SA1035 = list(set(df_SA1035U.clone_id.unique()).union(set(df_SA1035T.clone_id.unique())))
    # fill in missing clones with 0s
    df_SA1035T = fill_in_missing_clones(df_SA1035T, clone_list_SA1035)
    df_SA1035U = fill_in_missing_clones(df_SA1035U, clone_list_SA1035)

    # fill in the missing clones for SA535
    # find the set of clone_ids that appear in the union of the treated and untreated samples
    clone_list_SA535 = list(set(df_SA535U.clone_id.unique()).union(set(df_SA535T.clone_id.unique())))
    # fill in missing clones with 0s
    df_SA535T = fill_in_missing_clones(df_SA535T, clone_list_SA535)
    df_SA535U = fill_in_missing_clones(df_SA535U, clone_list_SA535)
    
    # sort timepoints again now that we've added the earliest untreated timepoint and filled in missing clones
    df_SA1035U = sort_timepoints(df_SA1035U)
    df_SA1035T = sort_timepoints(df_SA1035T)
    df_SA535U = sort_timepoints(df_SA535U)
    df_SA535T = sort_timepoints(df_SA535T)

    # compute the proxy s coefficients for all the samples
    proxy_s_df_SA535U = calc_proxy_s(df_SA535U)
    proxy_s_df_SA535T = calc_proxy_s(df_SA535T)
    proxy_s_df_SA1035U = calc_proxy_s(df_SA1035U)
    proxy_s_df_SA1035T = calc_proxy_s(df_SA1035T)

    # add columns to all proxy dfs to indicate the treatment status
    proxy_s_df_SA1035T['treatment'] = 'Rx+'
    proxy_s_df_SA1035U['treatment'] = 'Rx-'
    proxy_s_df_SA535T['treatment'] = 'Rx+'
    proxy_s_df_SA535U['treatment'] = 'Rx-'

    # combine the two dfs
    proxy_s_df_SA1035 = pd.concat([proxy_s_df_SA1035U, proxy_s_df_SA1035T], ignore_index=True)
    proxy_s_df_SA535 = pd.concat([proxy_s_df_SA535U, proxy_s_df_SA535T], ignore_index=True)

    # create a column for the relative ranking of proxy_s values within a treatment group
    # lowest ranking is the least fit clone, highest number is the most fit clone
    for treatment, chunk in proxy_s_df_SA1035.groupby('treatment'):
        proxy_s_df_SA1035.loc[proxy_s_df_SA1035.treatment==treatment, 'proxy_s_rank'] = chunk.proxy_s.rank(ascending=True)
    for treatment, chunk in proxy_s_df_SA535.groupby('treatment'):
        proxy_s_df_SA535.loc[proxy_s_df_SA535.treatment==treatment, 'proxy_s_rank'] = chunk.proxy_s.rank(ascending=True)
    
    # plot the proxy s values for all datasets
    plot_proxy_s_coefficients(proxy_s_df_SA1035, proxy_s_df_SA535, argv)


if __name__ == '__main__':
    main()
