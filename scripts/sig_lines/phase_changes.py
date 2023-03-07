import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_phase_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_g', type=str, help='full df for filtered G1/2-phase cells')
    p.add_argument('cn_s', type=str, help='full df for filtered S-phase cells')
    p.add_argument('cn_low', type=str, help='full df for filtered LQity cells')
    p.add_argument('cn_g_init', type=str, help='df for cells initialized as G1/2-phase')
    p.add_argument('cn_s_init', type=str, help='df for cells initialized as S-phase')
    p.add_argument('dataset', type=str, help='name of this dataset')
    p.add_argument('frac_rt_col', type=str, help='name of the column containing the fraction of replicated bins')
    p.add_argument('out_tsv', type=str, help='per-cell summary df for phase calls by all models and relevant features')
    p.add_argument('plot1', type=str, help='confusion matrix comparing the phase calls of all models')
    p.add_argument('plot2', type=str, help='violin plots of the features across different phase calls')

    return p.parse_args()


def plot_confusion_matrix(cn, argv):
    """ Plot a confusion matrix comparing laks_phase to pert_phase for each cell. """
    # create a figure
    fig, ax = plt.subplots(figsize=(4, 4))

    # plot a confusion matrix of cn where the rows are laks phase and the columns are pert phase
    sns.heatmap(pd.crosstab(cn['pert_phase'], cn['laks_phase']), annot=True, fmt='d', cmap='Blues', ax=ax)

    # rename the x and y axis labels
    ax.set_xlabel('Laks phase')
    ax.set_ylabel('PERT phase')
    ax.set_title(f'{argv.dataset}\n # cells')

    # save the figure
    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def plot_violin_plots(cn, argv):
    """ Use violin plots to show the distribution of features across different phase calls. """
    # list of columns to plot on the y-axis
    y_cols = [
        'corrected_breakpoints', 'corrected_madn', 'rpm_auto_norm',
        'rep_auto_norm', 'cell_frac_rep', 'quality',
        'gc_intercept', 'gc_slope', 'ploidy'
    ]

    # create a dictionary that maps each y-column name to a y-axis label
    y_label_dict = {
        'corrected_breakpoints': 'CN breakpoints',
        'corrected_madn': 'RPM median absolute deviation\nbetween neighboring bins',
        'rpm_auto_norm': 'RPM autocorrelation',
        'rep_auto_norm': 'PERT rep state autocorrelation',
        'cell_frac_rep': 'PERT fraction of replicated bins',
        'quality': 'Laks et al quality score',
        'gc_intercept': 'GC bias intercept',
        'gc_slope': 'GC bias slope',
        'ploidy': 'HMMcopy ploidy'
    }

    # create a matplotlib figure with 3 rows and 3 columns
    fig, axes = plt.subplots(3, 3, figsize=(12, 12), tight_layout=True)
    ax = axes.flatten()

    phase_cmap = get_phase_cmap()

    # loop through each column in y_cols, plotting the violin plot on the corresponding axis
    for y_col in y_cols:
        sns.violinplot(data=cn, x='laks_phase', hue='pert_phase', y=y_col, ax=ax[y_cols.index(y_col)], palette=phase_cmap, linewidth=1)
        ax[y_cols.index(y_col)].set_xlabel('Laks et al phase')
        ax[y_cols.index(y_col)].set_ylabel(y_label_dict[y_col])
        ax[y_cols.index(y_col)].legend(loc='upper right', title='PERT phase')
        ax[y_cols.index(y_col)].set_title(argv.dataset)

    # save the figure
    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    # load the filtered dataframes
    cn_g = pd.read_csv(argv.cn_g, sep='\t')
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_lowqual = pd.read_csv(argv.cn_low, sep='\t')

    # compute the ploidy and gc bias coefficients for each cell
    for cn in [cn_g, cn_s, cn_lowqual]:
        for cell_id, cell_cn in cn.groupby('cell_id'):
            # ploidy is the modal state
            cn.loc[cn.cell_id == cell_id, 'ploidy'] = cell_cn.state.value_counts().index[0]

            # fit a 1st order polynomial to rpm vs gc
            gc_fit = np.polyfit(cell_cn.gc, cell_cn.rpm, 1)
            cn.loc[cn.cell_id == cell_id, 'gc_slope'] = gc_fit[0]
            cn.loc[cn.cell_id == cell_id, 'gc_intercept'] = gc_fit[1]
    
    # subset to just cell-specific columns
    cell_columns = [
        'cell_id', 'total_mapped_reads_hmmcopy', 'is_s_phase', 'is_s_phase_prob', 'quality',
        'in_tree', 'clone_id', 'madn', 'lrs', 'corrected_madn', 'corrected_breakpoints', 'assigned_clone_id',
        'model_tau', 'model_u', 'cell_frac_rep', 'rpm_auto', 'rep_auto', 'cn_bk', 'rep_bk', 
        'frac_cn0', 'rpm_auto_norm', 'rep_auto_norm', 'gc_slope', 'gc_intercept', 'ploidy'
    ]

    cn_s = cn_s[cell_columns].drop_duplicates().reset_index(drop=True)
    cn_g = cn_g[cell_columns].drop_duplicates().reset_index(drop=True)
    cn_lowqual = cn_lowqual[cell_columns].drop_duplicates().reset_index(drop=True)

    # add column to both dataframes to indicate whether the cell is in G1/2 or S phase
    cn_g['pert_phase'] = 'G1/2'
    cn_s['pert_phase'] = 'S'
    cn_lowqual['pert_phase'] = 'LQ'

    # concatenate the two dataframes
    cn = pd.concat([cn_g, cn_s, cn_lowqual], ignore_index=True)

    # rename is_s_phase column to is_s_phase_laks
    cn = cn.rename(columns={'is_s_phase': 'is_s_phase_laks'})

    # load the dataframes for the cells initialized as G1/2-phase and S-phase
    cn_g_init = pd.read_csv(argv.cn_g_init, sep='\t')
    cn_s_init = pd.read_csv(argv.cn_s_init, sep='\t')

    # add columns to both dataframes to indicate whether the cell is in G1/2 or S phase
    cn_g_init['is_s_phase_init'] = 'G1/2'
    cn_s_init['is_s_phase_init'] = 'S'

    # subset to just the is_s_phase_init and cell_id columns
    cn_g_init = cn_g_init[['is_s_phase_init', 'cell_id']].drop_duplicates().reset_index(drop=True)
    cn_s_init = cn_s_init[['is_s_phase_init', 'cell_id']].drop_duplicates().reset_index(drop=True)

    # concatenate the two dataframes
    cn_init = pd.concat([cn_g_init, cn_s_init], ignore_index=True)

    # merge cn_init with cn
    cn = cn.merge(cn_init, on='cell_id', how='left')

    # create a new column named 'laks_phase' which has values 'G1/2', 'S', or 'LQ'
    # create a map of is_s_phase_laks to laks_phase
    cn['laks_phase'] = cn['is_s_phase_laks'].map({False: 'G1/2', True: 'S'})
    # for rows where quality < 0.75 and is_s_phase_laks is False, set laks_phase to 'LQ'
    cn.loc[(cn['quality'] < 0.75) & (cn['is_s_phase_laks'] == False), 'laks_phase'] = 'LQ'

    # plot the confusion matrix
    plot_confusion_matrix(cn, argv)

    # plot the violin plots
    plot_violin_plots(cn, argv)

    # save the dataframe
    cn.to_csv(argv.out_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
