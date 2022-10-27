import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_cell_cn_profile
from matplotlib.colors import ListedColormap
from scgenome import refgenome
from sklearn import preprocessing
from scipy.stats import hypergeom
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', type=str, help='full df for S-phase cells')
    p.add_argument('cn_g', type=str, help='full df for G1/2-phase cells')
    p.add_argument('cn_gr', type=str, help='full df for G1/2-phase cells that were passed to PERT as S-phase and subsequently recovered')
    p.add_argument('rt', type=str, help='table of RT pseudobulk profiles')
    p.add_argument('rep_col', type=str, help='column for replication state (relevant for RT pseudobulk profile column names)')
    p.add_argument('dataset', type=str, help='name of this dataset')
    p.add_argument('out_tsv', type=str, help='Table of the number of cells per cell cycle phase and clone')
    p.add_argument('plot1', type=str, help='figure comparing the clone pseudobulk rt profiles')
    p.add_argument('plot2', type=str, help='figure showing the distribution of cell cycle phases for each clone')

    return p.parse_args()


def plot_cell_cn_profile2(ax, cn_data, value_field_name, cn_field_name=None, max_cn=13,
                          chromosome=None, s=5, squashy=False, color=None, alpha=1,
                          lines=False, label=None, scale_data=False):
    """ Plot copy number profile on a genome axis

    Args:
        ax: matplotlib axis
        cn_data: copy number table
        value_field_name: column in cn_data to use for the y axis value
    
    Kwargs:
        cn_field_name: state column to color scatter points
        max_cn: max copy number for y axis
        chromosome: single chromosome plot
        s: size of scatter points

    The cn_data table should have the following columns (in addition to value_field_name and
    optionally cn_field_name):
        - chr
        - start
        - end
    """
    chromosome_info = refgenome.info.chromosome_info[['chr', 'chromosome_start', 'chromosome_end']].copy()
    chromosome_info['chr'] = pd.Categorical(chromosome_info['chr'], categories=cn_data['chr'].cat.categories)
    plot_data = cn_data.merge(chromosome_info)
    plot_data = plot_data[plot_data['chr'].isin(refgenome.info.chromosomes)]
    plot_data['start'] = plot_data['start'] + plot_data['chromosome_start']
    plot_data['end'] = plot_data['end'] + plot_data['chromosome_start']

    squash_coeff = 0.15
    squash_f = lambda a: np.tanh(squash_coeff * a)
    if squashy:
        plot_data[value_field_name] = squash_f(plot_data[value_field_name])
    
    if scale_data:
        plot_data[value_field_name] = preprocessing.scale(plot_data[value_field_name].values)

    if lines:
        chr_order = [str(i+1) for i in range(22)]
        chr_order.append('X')
        chr_order.append('Y')
        plot_data.chr.cat.set_categories(chr_order, inplace=True)
        plot_data = plot_data.sort_values(by=['chr', 'start'])
        if cn_field_name is not None:
            ax.plot(
                plot_data['start'], plot_data[value_field_name], alpha=0.3, c='k', label=''
            )
        elif color is not None:
            ax.plot(
                plot_data['start'], plot_data[value_field_name], alpha=0.3, c=color, label=''
            )
        else:
            ax.plot(
                plot_data['start'], plot_data[value_field_name], alpha=0.3, label=''
            )
    
    if label is None:
        label = value_field_name
    
    if cn_field_name is not None:
        ax.scatter(
            plot_data['start'], plot_data[value_field_name],
            c=plot_data[cn_field_name], s=s, alpha=alpha, label=label,
            cmap=get_cn_cmap(plot_data[cn_field_name].astype(int).values),
        )
    elif color is not None:
         ax.scatter(
            plot_data['start'], plot_data[value_field_name],
            c=color, s=s, alpha=alpha, label=label
        )
    else:
        ax.scatter(
            plot_data['start'], plot_data[value_field_name], s=s, alpha=alpha, label=label
        )

    if chromosome is not None:
        chromosome_length = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_length']
        chromosome_start = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_start']
        chromosome_end = refgenome.info.chromosome_info.set_index('chr').loc[chromosome, 'chromosome_end']
        xticks = np.arange(0, chromosome_length, 2e7)
        xticklabels = ['{0:d}M'.format(int(x / 1e6)) for x in xticks]
        xminorticks = np.arange(0, chromosome_length, 1e6)
        ax.set_xlabel(f'chromosome {chromosome}')
        ax.set_xticks(xticks + chromosome_start)
        ax.set_xticklabels(xticklabels)
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(xminorticks + chromosome_start))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlim((chromosome_start, chromosome_end))

    else:
        ax.set_xlim((-0.5, refgenome.info.chromosome_end.max()))
        ax.set_xlabel('chromosome')
        ax.set_xticks([0] + list(refgenome.info.chromosome_end.values))
        ax.set_xticklabels([])
        ax.xaxis.tick_bottom()
        ax.yaxis.tick_left()
        ax.xaxis.set_minor_locator(matplotlib.ticker.FixedLocator(refgenome.info.chromosome_mid))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.FixedFormatter(refgenome.info.chromosomes))

    if squashy:
        yticks = np.array([0, 2, 4, 7, 20])
        yticks_squashed = squash_f(yticks)
        ytick_labels = [str(a) for a in yticks]
        ax.set_yticks(yticks_squashed)
        ax.set_yticklabels(ytick_labels)
        ax.set_ylim((-0.01, 1.01))
        ax.spines['left'].set_bounds(0, 1)
    elif max_cn is not None:
        ax.set_ylim((-0.05*max_cn, max_cn))
        ax.set_yticks(range(0, int(max_cn) + 1))
        ax.spines['left'].set_bounds(0, max_cn)
    

    if chromosome is not None:
        sns.despine(ax=ax, offset=10, trim=False)
    else:
        sns.despine(ax=ax, offset=10, trim=True)

    return chromosome_info


def plot_clone_rt_profiles(rt, argv):
    cols = [x for x in rt.columns if x.startswith('pseduobulk_clone') and x.endswith(argv.rep_col)]
    clones = [x.split('_')[1].replace('clone', '') for x in cols]

    fig, ax = plt.subplots(4, 1, figsize=(16,16), tight_layout=True)
    ax = ax.flatten()
    i = 0
    for clone_id, col in zip(clones, cols):
        # plot the whole genome
        plot_cell_cn_profile2(
            ax[0], rt, col, color='C{}'.format(i), 
            max_cn=None, scale_data=False, lines=True, label=clone_id
        )
        # zoom in on chr1
        plot_cell_cn_profile2(
            ax[1], rt, col, color='C{}'.format(i), chromosome='1',
            max_cn=None, scale_data=False, lines=True, label=clone_id
        )
        
        # compute distance between this clone's RT and the reference clone (A)
        if clone_id != 'A':
            rt['temp_diff_rt'] = rt[col] - rt['pseduobulk_cloneA_{}'.format(argv.rep_col)]

            # plot the whole genome
            plot_cell_cn_profile2(
                ax[2], rt, 'temp_diff_rt', color='C{}'.format(i), 
                max_cn=None, scale_data=False, lines=True, label=clone_id
            )
            # zoom in on chr1
            plot_cell_cn_profile2(
                ax[3], rt, 'temp_diff_rt', color='C{}'.format(i), chromosome='1',
                max_cn=None, scale_data=False, lines=True, label=clone_id
            )
            
        
        i += 1

    for i in range(4):
        ax[i].set_title(argv.dataset)
        ax[i].legend(title='Clone ID')
        if i < 2:
            ax[i].set_ylabel('RT profile\n<--late | early-->')
        else:
            ax[i].set_ylabel('Relative RT to clone A\n<--clone later | clone earlier-->')

    fig.savefig(argv.plot1, bbox_inches='tight')


def clone_spf_analysis(cn_s, cn_g, cn_gr, argv):
    # add column to denote phase of each df
    cn_s['cell_cycle'] = 'S'
    cn_gr['cell_cycle'] = 'G1/2 recovered'
    cn_g['cell_cycle'] = 'G1/2 tree'

    # rename the clone_id column to match the cells in the tree
    cn_s['clone_id'] = cn_s['assigned_clone_id']
    cn_gr['clone_id'] = cn_gr['assigned_clone_id']

    # concatenate all cells into one df but only store relevant columns
    coi = ['cell_id', 'cell_cycle', 'clone_id']
    df = pd.concat([cn_s[coi], cn_gr[coi], cn_g[coi]], ignore_index=True)
    df = df.sort_values(['clone_id', 'cell_cycle']).drop_duplicates(ignore_index=True)

    clones = df['clone_id'].unique()

    # count the number of cells belonging to each clone and cell cycle phase
    df2 = df[['cell_cycle', 'clone_id']].value_counts().to_frame(name='num_cells').reset_index().sort_values(['clone_id', 'cell_cycle']).reset_index(drop=True)
    
    # add all the counts for clones+phase combos with no cells
    absent_df = []
    for clone_id in clones:
        if clone_id not in df2.query("cell_cycle=='S'")['clone_id'].values:
            absent_df.append(pd.DataFrame({
                'cell_cycle': ['S'], 'clone_id': [clone_id], 'num_cells': [0]
            }))
        if clone_id not in df2.query("cell_cycle=='G1/2 recovered'")['clone_id'].values:
            absent_df.append(pd.DataFrame({
                'cell_cycle': ['G1/2 recovered'], 'clone_id': [clone_id], 'num_cells': [0]
            }))
        if clone_id not in df2.query("cell_cycle=='G1/2 tree'")['clone_id'].values:
            absent_df.append(pd.DataFrame({
                'cell_cycle': ['G1/2 tree'], 'clone_id': [clone_id], 'num_cells': [0]
            }))
    
    # concatenate into one dataframe
    if len(absent_df) > 0:
        absent_df = pd.concat(absent_df, ignore_index=True)
        df2 = pd.concat([df2, absent_df], ignore_index=True)

    # create arrays of clone counts within each cell cycle phase
    num_cells_s = df2.query("cell_cycle=='S'")['num_cells'].values
    num_cells_g_tree = df2.query("cell_cycle=='G1/2 tree'")['num_cells'].values
    num_cells_g_rec = df2.query("cell_cycle=='G1/2 recovered'")['num_cells'].values
    # treat both the tree and recovered G1/2-phase cells as being nonreplicating here
    num_cells_g = num_cells_g_tree + num_cells_g_rec

    # convert clone counts to clone frequencies within each cell cycle phase
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

    # divide by total number of hypotheses tested (Bonferroni correction)
    positive_pvals *= 2 * len(clones)
    negative_pvals *= 2 * len(clones)

    fig, ax = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)

    # draw scatterplot comparing the relative fraction of each clone in S vs G1/2 phases
    pthresh = 1e-2
    for i, clone_id in enumerate(clones):
        if positive_pvals[i] < pthresh:
            ax[0].scatter(x=clone_frac_g[i], y=clone_frac_s[i], c='C{}'.format(i), label=clone_id, marker='^')
        elif negative_pvals[i] < pthresh:
            ax[0].scatter(x=clone_frac_g[i], y=clone_frac_s[i], c='C{}'.format(i), label=clone_id, marker='v')
        else:
            ax[0].scatter(x=clone_frac_g[i], y=clone_frac_s[i], c='C{}'.format(i), label=clone_id)

    # draw y=x line where we expect "neutral" clones to lie
    lims = [
        np.min([ax[0].get_xlim(), ax[0].get_ylim()]),  # min of both axes
        np.max([ax[0].get_xlim(), ax[0].get_ylim()]),  # max of both axes
    ]
    ax[0].plot(lims, lims, 'k--', alpha=0.25, zorder=0, label='neutral')

    ax[0].legend(title='Clone ID')
    ax[0].set_xlabel('Fraction of G1/2 cells assigned to clone')
    ax[0].set_ylabel('Fraction of S cells assigned to clone')
    ax[0].set_title('{}\nRelative proliferation rate of clones'.format(argv.dataset))

    # pie chart to represent the total fraction of cell cycle phases across all clones
    n_cells_per_phase = [sum(num_cells_g_tree), sum(num_cells_g_rec), sum(num_cells_s)]
    phases = ['G1/2 tree', 'G1/2 recovered', 'S']
    colors = matplotlib.cm.get_cmap('Pastel2', len(phases)).colors
    ax[1].pie(n_cells_per_phase, labels=phases, autopct='%1.1f%%', shadow=True, colors=colors)
    ax[1].set_title('Sample proliferation rate - {}'.format(argv.dataset))

    # save spf figure
    fig.savefig(argv.plot2, bbox_inches='tight')

    # save table used to generate spf figure
    df2.to_csv(argv.out_tsv, sep='\t', index=False)


def main():
    argv = get_args()
    # load long-form dataframes from different cell cycle phases
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g = pd.read_csv(argv.cn_g, sep='\t')
    cn_gr = pd.read_csv(argv.cn_gr, sep='\t')

    # load pseudobulk RT profiles
    rt = pd.read_csv(argv.rt, sep='\t')

    cn_s.chr = cn_s.chr.astype('str')
    cn_s.chr = cn_s.chr.astype('category')
    cn_g.chr = cn_g.chr.astype('str')
    cn_g.chr = cn_g.chr.astype('category')
    cn_gr.chr = cn_gr.chr.astype('str')
    cn_gr.chr = cn_gr.chr.astype('category')
    rt.chr = rt.chr.astype('str')
    rt.chr = rt.chr.astype('category')

    # merge end and gc columns into RT pseudobulks
    rt = pd.merge(rt, cn_s[['chr', 'start', 'end', 'gc']].drop_duplicates())

    # plot pseudoublk RT profiles
    plot_clone_rt_profiles(rt, argv)

    # plot S-phase fractions at the clone and sample levels
    # save both the plot and table used to make the plot
    clone_spf_analysis(cn_s, cn_g, cn_gr, argv)


if __name__ == '__main__':
    main()

