from argparse import ArgumentParser
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('-rr', '--rt_ref', type=str, help='rt pseudobulks when no labels are permuted')
    p.add_argument('-rp', '--rt_perm', type=str, nargs='+', help='rt pseudobulks for all the permuted datasets')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-rt', '--rt_table', help='table containing all the rt pseudobulk profiles for all permuted datasets')
    p.add_argument('-cp', '--cor_plot', help='plot showing the pairwise correlation between all permuted dataset bulk RT profiles')

    return p.parse_args()


def load_data(argv, i, d):
    rt = pd.read_csv(argv.rt_perm[i], sep='\t')
    rt['dataset'] = d
    return rt


def make_plots(legend_df, metrics_df, argv):
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    # barplot of the fraction of G1/2-cells accurately removed out of all those with swapped flow labels
    sns.barplot(data=legend_df, x='rate', y='accuracy', ax=ax[0])
    ax[0].set_ylabel('Fraction of mislabeled\ncells detected by model')
    ax[0].set_xlabel('Fraction of G1/2-phase cells mislabeled')
    ax[0].set_title('Model accuracy')

    # distribution of cell_frac_rep values based on the true flow sorting states
    sns.histplot(data=metrics_df.query("permuted==True"), x='cell_frac_rep', hue='true_cell_cycle_state', bins=20, multiple='stack')
    ax[1].set_title('Cells mislabeled as S-phase')
    ax[1].set_xlabel('Inferred fraction of replicated bins')

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

    # load rt bulks from dataset with no permuted labels (use as a reference)
    ref_rt = pd.read_csv(argv.rt_ref, sep='\t')

    # build table that matches config file for each permuted dataset
    legend_df = pd.DataFrame({
        'dataset': argv.datasets,
    })

    # load in all the permuted dataset bulk rt profiles
    all_rts = []
    i = 0
    for dataset, row in legend_df.groupby('dataset'):
        temp_rt = load_data(argv, i, dataset)
        all_rts.append(temp_rt)
        i += 1
    all_rts = pd.concat(all_rts, ignore_index=True)

    # merge each permuted dataset's relevant RT columns into rt_wide
    rt_wide = ref_rt.copy()
    bulk_rt_cols = ['rt_merged_T47D', 'rt_merged_GM18507']

    for dataset, temp_rt in all_rts.groupby('dataset'):
        # rename columns in temp_rt to reflect the current dataset
        t_col = 'T47D_{}'.format(dataset)
        gm_col = 'GM18507_{}'.format(dataset)
        bulk_rt_cols.extend([t_col, gm_col])
        temp_rt = temp_rt.rename(columns={
            'rt_T47D': t_col,
            'rt_GM18507': gm_col,
        })
        
        # merge the columns using loci
        temp_rt = temp_rt[['chr', 'start', t_col, gm_col]]
        rt_wide = pd.merge(rt_wide, temp_rt)
    
    # save the table used to compute the correlations
    rt_wide.to_csv(argv.rt_table, sep='\t', index=False)

    # remove columns that should not appear in the correlation matrix
    bad_columns = ['chr', 'start', 'end', 'gc', 'rt_split_T47D', 'rt_split_GM18507', 'rt_diff_split', 'rt_diff_merged']
    rt_wide = rt_wide.drop(columns=bad_columns)

    # rename the reference dataset columns with no permuted cells
    rt_wide = rt_wide.rename(columns={
        'rt_merged_T47D': 'T47D',
        'rt_merged_GM18507': 'GM18507',
    })

    # replace all the underscores with spaces in the column names
    rt_wide.columns = [col.replace('_', ' ') for col in rt_wide.columns]

    # compute the pairwise correlation between all bulk rt columns
    # sort the columns before computing the correlaiton to ensure the same order
    corr = rt_wide[sorted(rt_wide.columns)].corr()

    # plot the correlation matrix as a heatmap
    fig = plt.figure(figsize=(10, 10))

    # mask the upper triangle of the heatmap
    mask = np.zeros_like(corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    # plot the heatmap, including the mask and whitespace between the cells but not digit annotations
    sns.heatmap(corr, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=False)
    plt.title('Pseudobulk RT correlations\nPermuted datasets')

    # save the figure
    fig.savefig(argv.cor_plot, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
