import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('s_phase_rx', help='Cell counts by clone & time for S-phase Rx')
    p.add_argument('non_s_phase_rx', help='Cell counts by clone & time for G1-phase Rx')
    p.add_argument('s_phase_unrx', help='Cell counts by clone & time for S-phase UnRx')
    p.add_argument('non_s_phase_unrx', help='Cell counts by clone & time for G1-phase UnRx')
    p.add_argument('s_coeffs', help='Table of clone selection coefficients for all datasets')
    p.add_argument('dataset')
    p.add_argument('rx_dataset')
    p.add_argument('unrx_dataset')
    p.add_argument('output_pdf', help='output pdf comparing selection coefficients to SPF at the first unrx timepoint')
    p.add_argument('output_tsv', help='output df used to create plot')

    return p.parse_args()


def main():
    argv = get_args()

    s_rx_df = pd.read_csv(argv.s_phase_rx, sep='\t')
    g_rx_df = pd.read_csv(argv.non_s_phase_rx, sep='\t')
    s_unrx_df = pd.read_csv(argv.s_phase_unrx, sep='\t')
    g_unrx_df = pd.read_csv(argv.non_s_phase_unrx, sep='\t')

    # combine four input cell count dfs into one
    s_rx_df['is_s_phase'] = True
    s_rx_df['datasetname'] = argv.rx_dataset
    s_rx_df['Rx'] = False
    s_unrx_df['is_s_phase'] = True
    s_unrx_df['datasetname'] = argv.unrx_dataset
    s_unrx_df['Rx'] = False
    g_rx_df['is_s_phase'] = False
    g_rx_df['datasetname'] = argv.rx_dataset
    g_rx_df['Rx'] = True
    g_unrx_df['is_s_phase'] = False
    g_unrx_df['datasetname'] = argv.unrx_dataset
    g_unrx_df['Rx'] = False
    df = pd.concat([s_rx_df, s_unrx_df, g_rx_df, g_unrx_df], ignore_index=True)

    # load in fitClone selection coefficients and merge with df
    fitness_scores = pd.read_csv(argv.s_coeffs)
    fitness_scores.drop(columns=['timepoint', 'fraction', 'ncells', 'S'], inplace=True)
    fitness_scores.rename(columns={'clone_name': 'clone_id'}, inplace=True)
    fitness_scores.drop_duplicates(inplace=True)
    df = pd.merge(df, fitness_scores, on=['datasetname', 'clone_id'])

    # compute SPF for each clone using the UnRx data at the earliest timepoint
    first_time = min(df.timepoint.unique())
    first_unrx_df = df.loc[(df['timepoint']==first_time) & (df['Rx']==False)]
    early_unrx_spf_df = []
    for clone_id, chunk in first_unrx_df.groupby('clone_id'):
        s_row = chunk.query("is_s_phase==True")
        g_row = chunk.query("is_s_phase==False")
        num_s = s_row.num_cells.values[0]
        num_g = g_row.num_cells.values[0]
        spf = num_s / (num_s + num_g)
        t = chunk.timepoint.values[0]
        rx = chunk.Rx.values[0]
        label = chunk.labels.values[0]
        temp_df = pd.DataFrame({
            'timepoint': [t], 'clone_id': [clone_id], 'num_cells_s': [num_s], 'num_cells_g': [num_g], 
            'SPF': [spf], 'Rx': [rx], 'labels': [label]
        })
        early_unrx_spf_df.append(temp_df)
    early_unrx_spf_df = pd.concat(early_unrx_spf_df, ignore_index=True)

    # compute S-phase enrichment: the fraction of S-phase cells belonging to clone c
    # minus the fraction of G1/2-phase cells belonging to clone c at time t
    total_num_S = sum(early_unrx_spf_df['num_cells_s'].values)
    total_num_G = sum(early_unrx_spf_df['num_cells_g'].values)
    early_unrx_spf_df['SPE'] = (early_unrx_spf_df['num_cells_s'] / total_num_S) - (early_unrx_spf_df['num_cells_g'] / total_num_G)

    # merge early SPF with fitness selection coefficients
    s_coeffs = df[['clone_id', 'datasetname', 'labels', 'mean_1_plus_s', 'median_1_plus_s', 'sd_1_plus_s']].drop_duplicates()
    early_unrx_spf_df2 = early_unrx_spf_df[['clone_id', 'SPF', 'SPE', 'num_cells_s', 'num_cells_g']]
    early_unrx_spf_df3 = pd.merge(s_coeffs, early_unrx_spf_df2)
    early_unrx_spf_df3.reset_index(inplace=True, drop=True)

    # calculate the change in clone selection coefficients between UnRx and Rx data
    early_unrx_spf_df4 = []
    for clone_id, chunk in early_unrx_spf_df3.groupby('clone_id'):
        rx_row = chunk.loc[chunk['datasetname']==argv.rx_dataset]
        unrx_row = chunk.loc[chunk['datasetname']==argv.unrx_dataset]
        mean_1_plus_s_rx = rx_row['mean_1_plus_s'].values[0]
        median_1_plus_s_rx = rx_row['median_1_plus_s'].values[0]
        sd_1_plus_s_rx = rx_row['sd_1_plus_s'].values[0]
        mean_1_plus_s_unrx = unrx_row['mean_1_plus_s'].values[0]
        median_1_plus_s_unrx = unrx_row['median_1_plus_s'].values[0]
        sd_1_plus_s_unrx = unrx_row['sd_1_plus_s'].values[0]
        first_unrx_spf = unrx_row['SPF'].values[0]
        first_unrx_spe = unrx_row['SPE'].values[0]
        num_cells_g = unrx_row['num_cells_g'].values[0]

        temp_df = pd.DataFrame({
            'clone_id': [clone_id], 'first_unrx_spf': [first_unrx_spf], 'first_unrx_spe': [first_unrx_spe],
            'mean_1_plus_s_rx': [mean_1_plus_s_rx], 'median_1_plus_s_rx': [median_1_plus_s_rx],
            'sd_1_plus_s_rx': [sd_1_plus_s_rx], 'mean_1_plus_s_unrx': [mean_1_plus_s_unrx],
            'median_1_plus_s_unrx': [median_1_plus_s_unrx], 'sd_1_plus_s_unrx': [sd_1_plus_s_unrx],
            'num_cells_g': [num_cells_g]
        })
        early_unrx_spf_df4.append(temp_df)
    early_unrx_spf_df4 = pd.concat(early_unrx_spf_df4, ignore_index=True)
    early_unrx_spf_df4['unrx_minus_rx_median_1_plus_s'] = early_unrx_spf_df4['median_1_plus_s_unrx'] - early_unrx_spf_df4['median_1_plus_s_rx']
    early_unrx_spf_df4['unrx_minus_rx_mean_1_plus_s'] = early_unrx_spf_df4['mean_1_plus_s_unrx'] - early_unrx_spf_df4['mean_1_plus_s_rx']
    early_unrx_spf_df4['unrx_minus_rx_sd_1_plus_s'] = early_unrx_spf_df4['sd_1_plus_s_unrx'] + early_unrx_spf_df4['sd_1_plus_s_rx']

    # plot figures
    fig, ax = plt.subplots(2, 3, figsize=(12, 10), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(
        data=early_unrx_spf_df4, x='unrx_minus_rx_median_1_plus_s', y='first_unrx_spf', hue='clone_id', size='num_cells_g', ax=ax[0]
    )
    sns.scatterplot(
        data=early_unrx_spf_df4, x='median_1_plus_s_rx', y='first_unrx_spf', hue='clone_id', size='num_cells_g', ax=ax[1]
    )
    sns.scatterplot(
        data=early_unrx_spf_df4, x='median_1_plus_s_unrx', y='first_unrx_spf', hue='clone_id', size='num_cells_g', ax=ax[2]
    )
    ax[0].set_title(argv.dataset)
    ax[1].set_title(argv.rx_dataset)
    ax[2].set_title(argv.unrx_dataset)

    sns.scatterplot(
        data=early_unrx_spf_df4, x='unrx_minus_rx_median_1_plus_s', y='first_unrx_spe', hue='clone_id', size='num_cells_g', ax=ax[3]
    )
    sns.scatterplot(
        data=early_unrx_spf_df4, x='median_1_plus_s_rx', y='first_unrx_spe', hue='clone_id', size='num_cells_g', ax=ax[4]
    )
    sns.scatterplot(
        data=early_unrx_spf_df4, x='median_1_plus_s_unrx', y='first_unrx_spe', hue='clone_id', size='num_cells_g', ax=ax[5]
    )
    ax[3].set_title(argv.dataset)
    ax[4].set_title(argv.rx_dataset)
    ax[5].set_title(argv.unrx_dataset)

    fig.savefig(argv.output_pdf, bbox_inches='tight')

    # save output tsv
    early_unrx_spf_df4.to_csv(argv.output_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
