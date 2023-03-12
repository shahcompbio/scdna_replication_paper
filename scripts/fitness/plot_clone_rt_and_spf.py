import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib import cm
from scipy.stats import hypergeom
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.plot_utils import plot_cell_cn_profile2
from common.colors import get_clone_cmap, get_rx_cmap



def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', type=str, help='full df for S-phase cells')
    p.add_argument('cn_g', type=str, help='full df for G1/2-phase cells')
    p.add_argument('rt', type=str, help='table of RT pseudobulk profiles')
    p.add_argument('rep_col', type=str, help='column for replication state (relevant for RT pseudobulk profile column names)')
    p.add_argument('dataset', type=str, help='name of this dataset')
    p.add_argument('out_tsv', type=str, help='Table of the number of cells per cell cycle phase and clone')
    p.add_argument('plot1', type=str, help='figure comparing the clone pseudobulk rt profiles')
    p.add_argument('plot2', type=str, help='figure showing the distribution of cell cycle phases for each clone')

    return p.parse_args()


def plot_clone_rt_profiles(rt, argv):
    cols = [x for x in rt.columns if x.startswith('pseudobulk_clone') and x.endswith(argv.rep_col)]
    clones = [x.split('_')[1].replace('clone', '') for x in cols]
    ref_clone = clones[0]

    clone_cmap = get_clone_cmap()

    fig, ax = plt.subplots(4, 1, figsize=(16,16), tight_layout=True)
    ax = ax.flatten()
    i = 0
    for clone_id, col in zip(clones, cols):
        # plot the whole genome
        plot_cell_cn_profile2(
            ax[0], rt, col, color=clone_cmap[clone_id], 
            max_cn=None, scale_data=False, lines=True, label=clone_id
        )
        # zoom in on chr1
        plot_cell_cn_profile2(
            ax[1], rt, col, color=clone_cmap[clone_id], chromosome='1',
            max_cn=None, scale_data=False, lines=True, label=clone_id
        )
        
        # compute distance between this clone's RT and the reference clone (A)
        if clone_id != ref_clone:
            rt['temp_diff_rt'] = rt[col] - rt['pseudobulk_clone{}_{}'.format(ref_clone, argv.rep_col)]

            # plot the whole genome
            plot_cell_cn_profile2(
                ax[2], rt, 'temp_diff_rt', color=clone_cmap[clone_id], 
                max_cn=None, scale_data=False, lines=True, label=clone_id
            )
            # zoom in on chr1
            plot_cell_cn_profile2(
                ax[3], rt, 'temp_diff_rt', color=clone_cmap[clone_id], chromosome='1',
                max_cn=None, scale_data=False, lines=True, label=clone_id
            )
            
        
        i += 1

    for i in range(4):
        dataset = argv.dataset.replace('_CISPLATIN_Combined', '')
        ax[i].set_title(dataset)
        ax[i].legend(title='Clone ID')
        if i < 2:
            ax[i].set_ylabel('RT profile\n<--late | early-->')
        else:
            ax[i].set_ylabel('Relative RT to clone {}\n<--clone later | clone earlier-->'.format(ref_clone))

    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)



def compute_fracs_and_pvals(df):
    """
    For a given dataset, compute the S-phase fraction of each cell cycle phase and test
    each clone for enrichment or depletion of S-phase cells
    """
    clones = df['clone_id'].unique()
    
    num_cells_s = np.zeros(len(clones))
    num_cells_g = np.zeros(len(clones))
    for i, clone_id in enumerate(clones):
        num_cells_s[i] = df.query('cell_cycle=="S"').query('clone_id=="{}"'.format(clone_id))['num_cells'].values[0]
        num_cells_g[i] = df.query('cell_cycle=="G1/2"').query('clone_id=="{}"'.format(clone_id))['num_cells'].values[0]
        
    # convert clone counts to clone frequencies within each cell cycle phase
    print('num_cells_s', num_cells_s)
    print('num_cells_g', num_cells_g)
    clone_frac_s = num_cells_s / sum(num_cells_s)
    clone_frac_g = num_cells_g / sum(num_cells_g)
    
    # statistical test to see which clones are enriched/depleted for S-phase cells
    positive_pvals = np.zeros(len(clones))
    for i, clone_id in enumerate(clones):

        x = num_cells_s[i] # number of S-phase cells belonging to this clone
        m = sum(num_cells_s) # total number of S-phase cells
        n = sum(num_cells_g) # total number of G1/2-phase cells
        k = int(clone_frac_g[i] * (n + m)) # expected number of G1/2 + S phase cells belonging to this clone
        N = m + n  # total number of cells in entire population
        # use hypergeometric survival function to see if this clone has
        # more S-phase cells than expected (positively selected)
        positive_pvals[i] = hypergeom(M=N, n=m, N=k).sf(x)  

    # subtract positive pval from 1 to see if clone has 
    # significantly fewer S-phase cells than expected
    negative_pvals = 1 - positive_pvals
    
    # create a dataframe with one entry per clone with all relevant stats
    df_out = pd.DataFrame({
        'clone_id': clones,
        'num_cells_s': num_cells_s,
        'num_cells_g': num_cells_g,
        'clone_frac_s': clone_frac_s,
        'clone_frac_g': clone_frac_g,
        'positive_p': positive_pvals,
        'negative_p': negative_pvals
    })
    
    return df_out


def timepoint_wrapper_fracs_and_pvals(df):
    ''' Compute clone cell cycle fractions within each timepoint '''
    df_out = []
    for timepoint, chunk in df.groupby('timepoint'):
        print('timepoint', timepoint)
        print('chunk', chunk, sep='\n')
        temp_out = compute_fracs_and_pvals(chunk)
        temp_out['timepoint'] = timepoint
        df_out.append(temp_out)
    df_out = pd.concat(df_out, ignore_index=True)
    return df_out


def clone_spf_analysis(cn_s, cn_g, argv):
    # add column to denote phase of each df
    cn_s['cell_cycle'] = 'S'
    cn_g['cell_cycle'] = 'G1/2'

    # # rename the clone_id column to match the cells in the tree
    # cn_s['clone_id'] = cn_s['assigned_clone_id']

    # concatenate all cells into one df but only store relevant columns
    coi = ['cell_id', 'cell_cycle', 'clone_id', 'timepoint']
    df = pd.concat([cn_s[coi], cn_g[coi]], ignore_index=True)
    df = df.sort_values(['clone_id', 'cell_cycle', 'timepoint']).drop_duplicates(ignore_index=True)

    # count the number of cells belonging to each clone and cell cycle phase
    df2 = df[['cell_cycle', 'clone_id', 'timepoint']].value_counts().to_frame(name='num_cells').reset_index().sort_values(['clone_id', 'cell_cycle']).reset_index(drop=True)
    
    # remove clones with "None" ID
    df2 = df2.loc[df2['clone_id']!='None']

    clones = df2['clone_id'].unique()
    timepoints = df2['timepoint'].unique()

    # add all the counts for clones+phase combos with no cells
    absent_df = []
    for timepoint in timepoints:
        for clone_id in clones:
            if clone_id not in df2.query("cell_cycle=='S'").query("timepoint=='{}'".format(timepoint))['clone_id'].values:
                absent_df.append(pd.DataFrame({
                    'cell_cycle': ['S'], 'clone_id': [clone_id], 'timepoint': [timepoint], 'num_cells': [0]
                }))
            if clone_id not in df2.query("cell_cycle=='G1/2'").query("timepoint=='{}'".format(timepoint))['clone_id'].values:
                absent_df.append(pd.DataFrame({
                    'cell_cycle': ['G1/2'], 'clone_id': [clone_id], 'timepoint': [timepoint], 'num_cells': [0]
                }))
    
    # concatenate into one dataframe
    if len(absent_df) > 0:
        absent_df = pd.concat(absent_df, ignore_index=True)
        df2 = pd.concat([df2, absent_df], ignore_index=True)

    # compute cell cycle fractions per clone & timepoint along with p-values
    print('DataFrame prior to calculating p-values', df2, sep='\n')
    df2 = timepoint_wrapper_fracs_and_pvals(df2)

    # Bonferroni correction of p-values by the total number of hypotheses tested
    # x2 because we're doing a two-sided t-test
    df2['positive_p_adj'] = df2['positive_p'] * 2 * df2.shape[0]
    df2['negative_p_adj'] = df2['negative_p'] * 2 * df2.shape[0]

    # create column named timepoint_int that strips the 'X' prefix from the timepoint column and converts to integer
    df2['timepoint_int'] = df2['timepoint'].str.replace('X', '').astype(int)
    max_time = df2['timepoint_int'].max()
    min_time = df2['timepoint_int'].min()

    fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True)

    pthresh = 1e-2

    # create custom legend for clones & timepoints
    clone_cmap = get_clone_cmap()
    viridis = cm.get_cmap('viridis', 256)
    timepoint_cmap = {}
    clone_legend_elements = [
        Line2D([0], [0], marker='^', color='w', label='enriched', markerfacecolor='k', markersize=10),
        Line2D([0], [0], marker='v', color='w', label='depleted', markerfacecolor='k', markersize=10)
    ]
    timepoint_legend_elements = clone_legend_elements.copy()
    for i, c in enumerate(sorted(df2.clone_id.unique())):
        color = clone_cmap[c]
        clone_legend_elements.append(Patch(facecolor=color, label=c))

    for i, t in enumerate(sorted(df2.timepoint_int.unique())):
        color = viridis((t - min_time) / (max_time - min_time))
        timepoint_str = 'X{}'.format(t)
        timepoint_cmap[timepoint_str] = color
        timepoint_legend_elements.append(Patch(facecolor=color, label=t))

    # draw scatterplot comparing the relative fraction of each clone in S vs G1/2 phases
    for i, row in df2.iterrows():
        clone_id = row['clone_id']
        timepoint = row['timepoint']
        if row['positive_p_adj'] < pthresh:
            ax[0].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id], marker='^')
            ax[1].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=timepoint_cmap[timepoint], marker='^')
        elif row['negative_p_adj'] < pthresh:
            ax[0].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id], marker='v')
            ax[1].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=timepoint_cmap[timepoint], marker='v')
        else:
            ax[0].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id])
            ax[1].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=timepoint_cmap[timepoint])

    # draw y=x line where we expect "neutral" clones to lie
    lims = [
        np.min([ax[0].get_xlim(), ax[0].get_ylim()]),  # min of both axes
        np.max([ax[0].get_xlim(), ax[0].get_ylim()]),  # max of both axes
    ]
    ax[0].plot(lims, lims, 'k--', alpha=0.25, zorder=0)
    ax[1].plot(lims, lims, 'k--', alpha=0.25, zorder=0)

    dataset = argv.dataset.replace('_CISPLATIN_Combined', '')

    ax[0].legend(handles=clone_legend_elements, title='Clone ID')
    ax[0].set_xlabel('Fraction of G1/2 cells in timepoint assigned to clone')
    ax[0].set_ylabel('Fraction of S cells in timepoint assigned to clone')
    ax[0].set_title('Phase enrichment of {} clones'.format(dataset))
    
    ax[1].legend(handles=timepoint_legend_elements, title='timepoint')
    ax[1].set_xlabel('Fraction of G1/2 cells in timepoint assigned to clone')
    ax[1].set_ylabel('Fraction of S cells in timepoint assigned to clone')
    ax[1].set_title('Phase enrichment of {} clones'.format(dataset))

    # save spf figure
    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)

    # save table used to generate spf figure
    df2.to_csv(argv.out_tsv, sep='\t', index=False)


def main():
    argv = get_args()
    # load long-form dataframes from different cell cycle phases
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g = pd.read_csv(argv.cn_g, sep='\t')

    # load pseudobulk RT profiles
    rt = pd.read_csv(argv.rt, sep='\t')

    cn_s.chr = cn_s.chr.astype('str')
    cn_s.chr = cn_s.chr.astype('category')
    cn_g.chr = cn_g.chr.astype('str')
    cn_g.chr = cn_g.chr.astype('category')
    rt.chr = rt.chr.astype('str')
    rt.chr = rt.chr.astype('category')

    # merge end and gc columns into RT pseudobulks
    rt = pd.merge(rt, cn_s[['chr', 'start', 'end', 'gc']].drop_duplicates())

    # plot pseudoublk RT profiles
    plot_clone_rt_profiles(rt, argv)

    # if the dataset starts with SA1035, rename clone F to A as they are nearly identical to one another
    if argv.dataset.startswith('SA1035'):
        cn_s.loc[cn_s.clone_id == 'F', 'clone_id'] = 'A'
        cn_g.loc[cn_g.clone_id == 'F', 'clone_id'] = 'A'

    # plot S-phase fractions at the clone and sample levels
    # save both the plot and table used to make the plot
    clone_spf_analysis(cn_s, cn_g, argv)


if __name__ == '__main__':
    main()

