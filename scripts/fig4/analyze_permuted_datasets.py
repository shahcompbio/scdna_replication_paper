from argparse import ArgumentParser
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('-cg', '--cn_good', type=str, nargs='+', help='long-form dataframes for all cells thought to be replicating by the pyro model')
    p.add_argument('-cb', '--cn_bad', type=str, nargs='+', help='long-form dataframes for all cells thought to be non-replicating by the pyro model')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-r', '--rates', type=float, nargs='+', help='permutation rate for each dataset')
    p.add_argument('-so', '--summary_output', help='table for the accruacy of each permuted dataset')
    p.add_argument('-mo', '--metrics_output', help='table containing per-cell metrics for cells in all permuted datasets')
    p.add_argument('-sp', '--summary_plots', help='plot showing the model accuracy for catching mislabeled cells')
    p.add_argument('-cp', '--ccc_plots', help='plot showing the ccc features for mislabeled cells')

    return p.parse_args()


def load_data(argv, i, d):
    cn_good = pd.read_csv(argv.cn_good[i], sep='\t')
    cn_good['filtered'] = False

    cn_bad = pd.read_csv(argv.cn_bad[i], sep='\t')
    cn_bad['filtered'] = True

    cn = pd.concat([cn_good, cn_bad], ignore_index=True)
    cn['dataset'] = d
    return cn


def make_plots(legend_df, metrics_df, argv):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    # barplot of the fraction of G1/2-cells accurately removed out of all those with swapped flow labels
    sns.barplot(data=legend_df, x='rate', y='accuracy', ax=ax[1])
    ax[1].set_ylabel('Fraction of mislabeled\ncells detected by model')
    ax[1].set_xlabel('Fraction of G1/2-phase cells mislabeled')
    ax[1].set_title('Model accuracy')

    # distribution of cell_frac_rep values based on the true flow sorting states
    sns.histplot(data=metrics_df.query("permuted==True"), x='cell_frac_rep', hue='true_cell_cycle_state', bins=20, multiple='stack', ax=ax[0])
    ax[0].set_title('Cells mislabeled as S-phase')
    ax[0].set_xlabel('Inferred fraction of replicated bins')

    fig.savefig(argv.summary_plots, bbox_inches='tight')

    # scatterplots of ccc features for the the mislabled cells, colored by whether the model thinks the cell is replicating (extreme_cell_frac)
    fig, ax = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(data=metrics_df.query("permuted==True"), x='is_s_phase_prob', y='quality', hue='extreme_cell_frac', alpha=0.5, ax=ax[0])
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='is_s_phase_prob', y='corrected_madn', hue='extreme_cell_frac', alpha=0.5, ax=ax[1])
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='is_s_phase_prob', y='corrected_breakpoints', hue='extreme_cell_frac', alpha=0.5, ax=ax[2])
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='corrected_madn', y='quality', hue='extreme_cell_frac', alpha=0.5, ax=ax[3])
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='corrected_breakpoints', y='quality', hue='extreme_cell_frac', alpha=0.5, ax=ax[4])
    sns.scatterplot(data=metrics_df.query("permuted==True"), x='corrected_madn', y='corrected_breakpoints', hue='extreme_cell_frac', alpha=0.5, ax=ax[5])

    for i in range(6):
        ax[i].set_title('Cells mislabeled as S-phase')

    fig.savefig(argv.ccc_plots, bbox_inches='tight')


def main():
    argv = get_args()

    # build table that matches config file for each permuted dataset
    # rates = [0.01] * 3
    # rates.extend([0.03]*3)
    # rates.extend([0.05]*3)
    # rates.extend([0.1]*3)
    # rates.extend([0.2]*3)
    # rates.extend([0.3]*3)
    legend_df = pd.DataFrame({
        'dataset': argv.datasets,
        'rate': argv.rates
    })

    # create a dict of DataFrames, one entry for each dataset
    all_cns = {}
    all_metrics = {}

    i = 0
    for dataset, row in legend_df.groupby('dataset'):
        temp_cn = load_data(argv, i, dataset)
        
        # add column to denote which cells have had their labels swapped for this dataset
        temp_cn['permuted'] = False
        temp_cn.loc[temp_cn['true_cell_cycle_state']!=temp_cn['cell_cycle_state'], 'permuted'] = True

        temp_metrics = temp_cn[[
            'cell_id', 'cell_frac_rep', 'true_cell_cycle_state', 'cell_cycle_state',
            'total_mapped_reads_hmmcopy', 'breakpoints', 'is_s_phase_prob',
            'madn', 'lrs', 'corrected_madn', 'corrected_breakpoints', 'quality',
            'extreme_cell_frac', 'filtered', 'permuted', 'dataset'
        ]].drop_duplicates().reset_index(drop=True)

        # what percentage of the cells with incorrect labels are properly identified by the model
        num_permuted = temp_metrics.query("permuted==True").shape[0]
        num_model_extreme = temp_metrics.query("permuted==True").query("extreme_cell_frac==True").shape[0]
        legend_df.loc[row.index, 'num_permuted'] = num_permuted
        legend_df.loc[row.index, 'num_model_extreme'] = num_model_extreme
        legend_df.loc[row.index, 'accuracy'] = num_model_extreme / num_permuted

        # save metrics and cn for this dataset in case I want to access them later
        all_cns[dataset] = temp_cn
        all_metrics[dataset] = temp_metrics

        i += 1


    # make one large metrics_df for all cells in all datasets
    metrics_df = []
    for k, v in all_metrics.items():
        metrics_df.append(v)
    metrics_df = pd.concat(metrics_df, ignore_index=True)
    metrics_df.head()


    make_plots(legend_df, metrics_df, argv)

    legend_df.to_csv(argv.summary_output, sep='\t', index=False)

    metrics_df.to_csv(argv.metrics_output, sep='\t', index=False)


if __name__ == '__main__':
    main()
