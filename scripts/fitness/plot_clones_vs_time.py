import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('treated_input', type=str, help='Table of the number of cells per cell cycle phase and clone for the treated sample')
    p.add_argument('untreated_input', type=str, help='Table of the number of cells per cell cycle phase and clone for the untreated sample')
    p.add_argument('dataset')
    p.add_argument('output_clones', type=str, help='figure showing the evolution of clones over time for each phase')
    p.add_argument('output_total', type=str, help='figure showing total number of cells at each timepoint')

    return p.parse_args()


def compute_clone_fracs_vs_time(df):
    """
    Computes the fraction of cells in each clone for each timepoint.
    """
    clone_frac_g_vs_time = np.zeros((len(df.clone_id.unique()), len(df.timepoint_int.unique())))
    clone_frac_s_vs_time = np.zeros((len(df.clone_id.unique()), len(df.timepoint_int.unique())))
    clone_legend = []
    timepoints = df.timepoint_int.unique()
    i = 0
    for clone_id, chunk in df.groupby('clone_id'):
        chunk.sort_values(by='timepoint_int', inplace=True)
        clone_frac_g_vs_time[i] = chunk['clone_frac_g'].values
        clone_frac_s_vs_time[i] = chunk['clone_frac_s'].values
        clone_legend.append(clone_id)
        i += 1
    return timepoints, clone_frac_g_vs_time, clone_frac_s_vs_time, clone_legend


def plot_clone_stackplots(df_U, df_T, argv):
    """
    Plots the evolution of clones over time for each cell cycle phase. 
    The top row shows the untreated sample, the bottom row shows the treated sample.
    """
    fig, ax = plt.subplots(2, 2, figsize=(8, 8), tight_layout=True, sharey=True)
    ax = ax.flatten()

    # untreated sample
    timepoints, clone_frac_g_vs_time, clone_frac_s_vs_time, clone_legend = compute_clone_fracs_vs_time(df_U)
    ax[0].stackplot(timepoints, clone_frac_g_vs_time, labels=clone_legend)
    ax[1].stackplot(timepoints, clone_frac_s_vs_time, labels=clone_legend)
    ax[0].set_ylabel('Clone fraction')
    ax[0].set_xlabel('Timepoint')
    ax[1].set_xlabel('Timepoint')
    ax[0].set_title('{}U: G1/2-phase'.format(argv.dataset))
    ax[1].set_title('{}U: S-phase'.format(argv.dataset))
    ax[1].legend(title='Clone ID')

    # treated sample
    timepoints, clone_frac_g_vs_time, clone_frac_s_vs_time, clone_legend = compute_clone_fracs_vs_time(df_T)
    ax[2].stackplot(timepoints, clone_frac_g_vs_time, labels=clone_legend)
    ax[3].stackplot(timepoints, clone_frac_s_vs_time, labels=clone_legend)
    ax[2].set_ylabel('Clone fraction')
    ax[2].set_xlabel('Timepoint')
    ax[3].set_xlabel('Timepoint')
    ax[2].set_title('{}T: G1/2-phase'.format(argv.dataset))
    ax[3].set_title('{}T: S-phase'.format(argv.dataset))
    ax[3].legend(title='Clone ID')

    fig.savefig(argv.output_clones, bbox_inches='tight', dpi=300)


def compute_total_cell_count_vs_time(df):
    """
    Computes the total number of cells in each cell cycle phase for each timepoint.
    """
    timepoints = df.timepoint_int.unique()
    total_counts = np.zeros((2, len(timepoints)))
    phase_legend = ['S', 'G1/2']
    i = 0
    for time, chunk in df.groupby('timepoint_int'):
        total_counts[0, i] = sum(chunk['num_cells_s'].values)
        total_counts[1, i] = sum(chunk['num_cells_g'].values)
        i += 1
    return timepoints, total_counts, phase_legend


def plot_total_cell_count_stackplot(df_U, df_T, argv):
    """
    Plots the total number of cells in each cell cycle phase at each timepoint.
    The top row shows the untreated sample, the bottom row shows the treated sample.
    """
    fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True)
    ax = ax.flatten()

    # untreated sample
    timepoints, total_counts, phase_legend = compute_total_cell_count_vs_time(df_U)
    ax[0].stackplot(timepoints, total_counts, labels=phase_legend)
    ax[0].set_ylabel('# cells')
    ax[0].set_xlabel('Timepoint')
    ax[0].set_title('{}U: total cell count'.format(argv.dataset))
    ax[0].legend(title='Phase')

    # treated sample
    timepoints, total_counts, phase_legend = compute_total_cell_count_vs_time(df_T)
    ax[1].stackplot(timepoints, total_counts, labels=phase_legend)
    ax[1].set_ylabel('# cells')
    ax[1].set_xlabel('Timepoint')
    ax[1].set_title('{}T: total cell count'.format(argv.dataset))
    ax[1].legend(title='Phase')

    fig.savefig(argv.output_total, bbox_inches='tight', dpi=300)


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


def main():
    argv = get_args()
    # load long-form dataframes from different cell cycle phases
    df_T = pd.read_csv(argv.treated_input, sep='\t')
    df_U = pd.read_csv(argv.untreated_input, sep='\t')

    # sort timepoints based on the timepoint_int column
    df_T = sort_timepoints(df_T)
    df_U = sort_timepoints(df_U)

    # find the earliest timepoint in the untreated sample
    earliest_timepoint = sorted(df_U.timepoint_int.unique())[0]
    # add data from the earliest untreated timepiont to the treated dataframe
    df_T = pd.concat([df_U[df_U.timepoint_int == earliest_timepoint], df_T], ignore_index=True)

    # find the set of clone_ids that appear in the union of the treated and untreated samples
    clone_list = list(set(df_U.clone_id.unique()).union(set(df_T.clone_id.unique())))
    # fill in missing clones with 0s
    df_T = fill_in_missing_clones(df_T, clone_list)
    df_U = fill_in_missing_clones(df_U, clone_list)

    # sort timepoints again now that we've added the earliest untreated timepoint and filled in missing clones
    df_T = sort_timepoints(df_T)
    df_U = sort_timepoints(df_U)

    # plot clone fractions for each phase & timepoint in the form of a stackplot
    plot_clone_stackplots(df_U, df_T, argv)

    # plot the total number of cells in each cell cycle phase for each timepoint
    # this should be helpful for identifying timepoints that should be removed
    plot_total_cell_count_stackplot(df_U, df_T, argv)


if __name__ == '__main__':
    main()
