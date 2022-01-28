import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('s_phase', help='input cell_id to clone_id tsv file for s phase cells')
    p.add_argument('non_s_phase', help='input cell_id to clone_id tsv file for non s phase cells')
    p.add_argument('times', help='time between passage numbers from fitness_time_scales.tsv')
    p.add_argument('datatag')
    p.add_argument('s_out', help='output of S-phase cell counts by time and clone')
    p.add_argument('non_s_out', help='output of non-S-phase cell counts by time and clone')
    p.add_argument('output_pdf', help='output pdf file fraction of each clone for S and non-S cells')

    return p.parse_args()


def main():
    argv = get_args()

    s_df = pd.read_csv(argv.s_phase, sep='\t')
    non_s_df = pd.read_csv(argv.non_s_phase, sep='\t')

    plot_label = s_df.datasetname.unique()[0] + ': ' + s_df.label.unique()[0]

    # load in timepoint to time (days) mapping for this dataset
    time_df = pd.read_csv(argv.times, sep='\t')
    time_df = time_df.loc[time_df['datatag']==argv.datatag]
    time_df = time_df[['timepoint', 'time']]


    # reduce to just unique mappings of cells to clones
    s_df = s_df[['cell_id', 'clone_id', 'timepoint']].drop_duplicates()
    non_s_df = non_s_df[['cell_id', 'clone_id', 'timepoint']].drop_duplicates()

    # find all timepoints & clones associated with this dataset
    times = sorted(non_s_df.timepoint.unique())
    times = [int(t.replace('X', '')) for t in times]
    clones = sorted(non_s_df.clone_id.unique())

    all_s_counts = []
    all_non_s_counts = []

    # find number of cells in each clone at each time
    for t in times:
        time = 'X{}'.format(t)
        s_chunk = s_df[s_df['timepoint'] == time]
        non_s_chunk = non_s_df[non_s_df['timepoint'] == time]

        # count nubmer of cells in each clone
        s_counts = s_chunk.clone_id.value_counts()
        non_s_counts = non_s_chunk.clone_id.value_counts()

        # sort both dicts alphabetically
        s_counts = dict(sorted(s_counts.items()))
        non_s_counts = dict(sorted(non_s_counts.items()))

        total_s_cells = sum(s_counts.values()) + np.finfo(float).eps
        total_non_s_cells = sum(non_s_counts.values()) + np.finfo(float).eps

        # loop through each clone and add number of cells to all counts dfs
        for i in clones:
            s_num_cells = s_counts[i] if i in s_counts else 0
            s_frac = s_num_cells / total_s_cells
            temp_s = pd.DataFrame({'timepoint': [t], 'clone_id': [i], 'num_cells': [s_num_cells], 'fraction': [s_frac]})
            non_s_num_cells = non_s_counts[i] if i in non_s_counts else 0
            non_s_frac = non_s_num_cells / total_non_s_cells
            temp_non_s = pd.DataFrame({'timepoint': [t], 'clone_id': [i], 'num_cells': [non_s_num_cells], 'fraction': [non_s_frac]})

            all_s_counts.append(temp_s)
            all_non_s_counts.append(temp_non_s)

    # contat list of dfs into one df
    all_s_counts = pd.concat(all_s_counts)
    all_non_s_counts = pd.concat(all_non_s_counts)

    # merge in real time (days) if time_df isn't empty
    if len(time_df) > 0:
        # remove X prefix and convert to int so it matches count dfs
        time_df['timepoint'] = time_df['timepoint'].apply(lambda x: int(x.replace('X', '')))
        # merge time_df into each count df
        all_non_s_counts = pd.merge(all_non_s_counts, time_df, on='timepoint')
        all_s_counts = pd.merge(all_s_counts, time_df, on='timepoint')
        # make sure the new time column is int64 just like other int columns
        time_col = 'time (days)'
        all_non_s_counts[time_col] = all_non_s_counts['time'].astype('int64')
        all_s_counts[time_col] = all_s_counts['time'].astype('int64')
    else:
        time_col = 'timepoint'

    # convert columns from objects to ints and floats
    all_s_counts['timepoint'] = all_s_counts['timepoint'].astype('int64')
    all_s_counts['num_cells'] = all_s_counts['num_cells'].astype('int64')
    all_s_counts['fraction'] = all_s_counts['fraction'].astype('float')
    all_non_s_counts['timepoint'] = all_non_s_counts['timepoint'].astype('int64')
    all_non_s_counts['num_cells'] = all_non_s_counts['num_cells'].astype('int64')
    all_non_s_counts['fraction'] = all_non_s_counts['fraction'].astype('float')

    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    axs = axs.flatten()

    # x-axis is time, y-axis is number of cells, hue is clone_id
    # do this for both S and non-S groups
    sns.scatterplot(data=all_s_counts, x=time_col, y='num_cells', hue='clone_id', ax=axs[0])
    sns.scatterplot(data=all_non_s_counts, x=time_col, y='num_cells', hue='clone_id', ax=axs[1])
    axs[0].set_title('{}\nS-phase'.format(plot_label))
    axs[1].set_title('{}\nnon-S-phase'.format(plot_label))

    # x-axis is time, y-axis is fraction of cells (normalized to this timept), hue is clone_id
    # do this for both S and non-S groups
    sns.scatterplot(data=all_s_counts, x=time_col, y='fraction', hue='clone_id', ax=axs[2])
    sns.scatterplot(data=all_non_s_counts, x=time_col, y='fraction', hue='clone_id', ax=axs[3])
    axs[2].set_title('{}\nS-phase'.format(plot_label))
    axs[3].set_title('{}\nnon-S-phase'.format(plot_label))

    fig.savefig(argv.output_pdf, bbox_inches='tight')

    all_s_counts.to_csv(argv.s_out, sep='\t', index=False)
    all_non_s_counts.to_csv(argv.non_s_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
