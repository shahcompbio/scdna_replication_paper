import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from scdna_replication_tools.plot_utils import plot_cell_cn_profile2, get_htert_cmap, get_cna_cmap
from common.colors import get_bkpt_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('-ir', '--input_rt', type=str, nargs='+', help='pseduobulk rt profiles for each dataset')
    p.add_argument('-ic', '--input_cn', type=str, nargs='+', help='pseduobulk cn profiles for each dataset')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='list of dataset names')
    p.add_argument('--plot1', help='plots showing the RT distribution against CNA types and breakpoints')
    p.add_argument('--plot2', help='plots reference RT and pseudobulk CN profiles')

    return p.parse_args()


def load_clone_rt_data(datasets, rt_paths):
    # load dataset pseudobulk rt profiles
    rt = pd.DataFrame()
    for d, rt_path in zip(datasets, rt_paths):
        temp_rt = pd.read_csv(rt_path, sep='\t')
        # set chr and start as the index columns
        temp_rt = temp_rt.set_index(['chr', 'start'])
        # subset to only the columns containing 'clone' and 'model_rep_state' substrings
        temp_rt = temp_rt[[c for c in temp_rt.columns if 'clone' in c and 'model_rep_state' in c]]
        # add the dataset d as the prefix for all the column names
        temp_rt = temp_rt.add_prefix('{}_'.format(d))
        # remove the 'pseudobulk_' substring from the column names
        temp_rt = temp_rt.rename(columns={c: c.replace('pseudobulk_', '') for c in temp_rt.columns})
        # replace the 'model_rep_state' substring with 'pseudobulk_rt' in the column names
        temp_rt = temp_rt.rename(columns={c: c.replace('model_rep_state', 'pseudobulk_rt') for c in temp_rt.columns})
        
        # reset the index
        temp_rt = temp_rt.reset_index()
        if rt.empty:
            rt = temp_rt
        else:
            rt = pd.merge(rt, temp_rt)

    # set chr column to category
    rt.chr = rt.chr.astype('str').astype('category')
    rt['end'] = rt['start'] + 500000 - 1
            
    return rt
    

def load_clone_cn_data(datasets, cn_paths):
    # load dataset pseudobulk cn profiles
    cn = pd.DataFrame()

    for d, cn_path in zip(datasets, cn_paths):
        temp_cn = pd.read_csv(cn_path, sep='\t')

        # set chr and start as the index columns
        temp_cn = temp_cn.set_index(['chr', 'start'])
        # subset to only the columns containing 'clone' and 'model_rep_state' substrings
        temp_cn = temp_cn[[c for c in temp_cn.columns if 'clone' in c]]
        # add the dataset d as the prefix for all the column names
        temp_cn = temp_cn.add_prefix('{}_'.format(d))
        # replace the 'clone_' substring with 'clone' in the column names
        temp_cn = temp_cn.rename(columns={c: c.replace('clone_', 'clone') for c in temp_cn.columns})
        # add the '_pseudobulk_cn' suffix to the column names
        temp_cn = temp_cn.rename(columns={c: '{}_pseudobulk_cn'.format(c) for c in temp_cn.columns})

        # reset the index
        temp_cn = temp_cn.reset_index()
        
        if cn.empty:
            cn = temp_cn
        else:
            cn = pd.merge(cn, temp_cn)

    # set chr column to category
    cn.chr = cn.chr.astype('str').astype('category')
    cn['end'] = cn['start'] + 500000 - 1

    return cn


def make_bk_clones(cn, clones):
    # create a table of copy number breakpoints for sample pseudobulk data
    # copy over the chromosome, start, and end columns
    cn_breakpoints = cn[['chr', 'start', 'end']].copy()
    for c in clones:
        # create a column for the copy number breakpoints
        cn_breakpoints['{}_pseudobulk_cn_breakpoints'.format(c)] = cn['{}_pseudobulk_cn'.format(c)].diff().fillna(0).astype('int')
        # convert all nonzero values to 1
        cn_breakpoints['{}_pseudobulk_cn_breakpoints'.format(c)] = cn_breakpoints['{}_pseudobulk_cn_breakpoints'.format(c)].apply(lambda x: 1 if x!=0 else 0)

    return cn_breakpoints


def compute_relative_rt_and_cn_clones(cn, rt, clones):
    for c in clones:
        ref_rt_col = 'SA039_cloneA_pseudobulk_rt'
        ref_cn_col = 'SA039_cloneA_pseudobulk_cn'
        
        temp_rt_col = '{}_pseudobulk_rt'.format(c)
        temp_cn_col = '{}_pseudobulk_cn'.format(c)
        
        relative_rt_col = '{}_pseudobulk_relative_rt'.format(c)
        relative_cn_col = '{}_pseudobulk_relative_cn'.format(c)
        
        rt[relative_rt_col] = rt[temp_rt_col] - rt[ref_rt_col]
        cn[relative_cn_col] = (cn[temp_cn_col] / cn[temp_cn_col].mode().values[0]) - (cn[ref_cn_col] / cn[ref_cn_col].mode().values[0])
    
    return cn, rt


def merge_cn_and_rt_info_clones(cn, rt, bk, clones):
    # merge the cn and rt tables into one long-form dataframe
    df = []

    for c in clones:
        ref_rt_col = 'SA039_cloneA_pseudobulk_rt'
        ref_cn_col = 'SA039_cloneA_pseudobulk_cn'
        
        temp_rt_col = '{}_pseudobulk_rt'.format(c)
        temp_cn_col = '{}_pseudobulk_cn'.format(c)
        temp_bk_col = '{}_pseudobulk_cn_breakpoints'.format(c)
        
        relative_rt_col = '{}_pseudobulk_relative_rt'.format(c)
        relative_cn_col = '{}_pseudobulk_relative_cn'.format(c)
        
        temp_df = pd.DataFrame({
            'chr': cn['chr'], 'start': cn['start'], 'end': cn['end'], 'clone': [c]*cn.shape[0],
            'WT_pseudobulk_rt': rt[ref_rt_col], 'WT_pseudobulk_cn': cn[ref_cn_col],
            'pseudobulk_rt': rt[temp_rt_col], 'pseudobulk_cn': cn[temp_cn_col], 'pseudobulk_bk': bk[temp_bk_col],
            'relative_rt': rt[relative_rt_col], 'relative_cn': cn[relative_cn_col]
        })
        
        df.append(temp_df)

    df = pd.concat(df, ignore_index=True)

    return df


def violins_with_pvals(df, x, y, hue, ax, box_pairs, order=None, test='t-test_ind', text_format='star', loc='inside', verbose=0, palette=None):
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, order=order, palette=palette)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test, order=order,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_WT_rt_vs_cna_type(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' 
    Plot the distribution of SA039 (WT) RT values against the ensuing CNAs that emerge
    in the other hTERT cell lines. This data will show whether gains and losses preferentially
    emerge in early or late replicating regions.
    '''
    cna_cmap = get_cna_cmap()
    x = "cna_type"
    y = "WT_pseudobulk_rt"
    hue = None
    box_pairs = [
        ('loss', 'gain'),
        ('loss', 'unaltered'),
        ('unaltered', 'gain'),
    ]
    order = ['loss', 'unaltered', 'gain']
    violins_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose, palette=cna_cmap)
    return ax


def plot_WT_rt_vs_bk(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' 
    Plot the distribution of SA039 (WT) RT values against the ensuing breakpoints that emerge
    in the other hTERT cell lines. This data will show whether CN breakpoints preferentially
    emerge in early or late replicating regions.
    '''
    bk_cmap = get_bkpt_cmap()
    x = "breakpoint"
    y = "WT_pseudobulk_rt"
    hue = None
    box_pairs = [
        ('No', 'Yes'),
    ]
    order = ['No', 'Yes']
    violins_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose, palette=bk_cmap)
    return ax



def plot_rt_distributions(df, argv):
    fig, ax = plt.subplots(2, 2, figsize=(8,8), tight_layout=True)
    ax = ax.flatten()

    # violin plots of RT distributions split by CNA type
    plot_WT_rt_vs_cna_type(df, ax[0])
    ax[0].set_xlabel('Subclonal CN in hTERT mutant cell lines')
    ax[0].set_ylabel('Ancestral hTERT RT\n<--late | early-->')
    ax[0].set_title('Location of subclonal CNAs')

    # histogram of RT values split by CNA type
    df['CNA type'] = df['cna_type']
    sns.histplot(
        data=df, x='WT_pseudobulk_rt', hue='CNA type', 
        common_norm=False, stat='density', ax=ax[1],
        palette=get_cna_cmap()
    )
    ax[1].set_xlabel('Ancestral hTERT RT\n<--late | early-->')
    ax[1].set_title('Location of subclonal CNAs')

    # violin plots of RT distributions split by breakpoint presence
    plot_WT_rt_vs_bk(df, ax[2])
    ax[2].set_xlabel('Subclonal CN breakpoints in hTERT mutant cell lines')
    ax[2].set_ylabel('Ancestral hTERT RT\n<--late | early-->')
    ax[2].set_title('Location of subclonal CNA breakpoints')

    # histogram of RT values split by breakpoint presence
    sns.histplot(
        data=df, x='WT_pseudobulk_rt', hue='breakpoint', 
        common_norm=False, stat='density', ax=ax[3],
        palette=get_bkpt_cmap()
    )
    ax[3].set_xlabel('Ancestral hTERT RT\n<--late | early-->')
    ax[3].set_title('Location of subclonal CNA breakpoints')

    # save the figure
    fig.savefig(argv.plot1, dpi=300, bbox_inches='tight')


def plot_profiles(cn, rt, clones, argv):
    # plot the reference RT profile and the relative CN profiles for each dataset
    # get the reference pseudobulk rt column
    ref_rt_col = 'SA039_cloneA_pseudobulk_rt'
    chrom_labels_to_remove = ['17', '19', '21']

    fig, ax = plt.subplots(2, 1, figsize=(16,8), tight_layout=True)
    ax = ax.flatten()

    # plot pseudobulk rt values of the reference WT dataset
    plot_cell_cn_profile2(
        ax[0], rt, ref_rt_col,
        max_cn=None, scale_data=False, lines=True,
        rasterized=True, rawy=True, s=1,
        chrom_labels_to_remove=chrom_labels_to_remove
    )
    ax[0].set_title('Ancestral hTERT RT (SA039 clone A)')
    ax[0].set_ylabel('Pseudobulk RT\n<--late | early-->')

    cmap = get_htert_cmap()
    for i, c in enumerate(clones):
        d = c.split('_')[0]
        relative_cn_col = '{}_pseudobulk_relative_cn'.format(c)
        # plot pseudobulk cn values for this dataset
        if 'cloneA' in c:
            plot_cell_cn_profile2(
                ax[1], cn, relative_cn_col, color=cmap[d], rawy=True, s=1,
                max_cn=None, scale_data=False, lines=True, label=d, rasterized=True,
                chrom_labels_to_remove=chrom_labels_to_remove
            )
        else:
            plot_cell_cn_profile2(
                ax[1], cn, relative_cn_col, color=cmap[d], rawy=True, s=1,
                max_cn=None, scale_data=False, lines=True, label='', rasterized=True,
                chrom_labels_to_remove=chrom_labels_to_remove
            )
    ax[1].set_title('Clone pseudobulk CN')
    ax[1].set_ylabel('ploidy-normalized CN\n<--losses | gains-->')
    ax[1].legend(ncol=2)

    # manually set the y-ticks for ax[1] to range from -1 to 1, spaced by 0.2
    ax[1].set_yticks(np.arange(-1, 5, 0.5))
    ax[1].spines['left'].set_bounds(-1, 4.5)

    fig.savefig(argv.plot2, dpi=300, bbox_inches='tight')


def main():
    argv = get_args()

    # load the rt and cn data
    # load the rt and cn data
    rt_clone = load_clone_rt_data(argv.datasets, argv.input_rt)
    cn_clone = load_clone_cn_data(argv.datasets, argv.input_cn)

    # remove the site of p53 deletion on chr17
    chr17_start_thresh = 21000001
    cn_clone = cn_clone.loc[~((cn_clone['chr'] == '17') & (cn_clone['start'] <= chr17_start_thresh))]
    rt_clone = rt_clone.loc[~((rt_clone['chr'] == '17') & (rt_clone['start'] <= chr17_start_thresh))]

    htert_clones = [c.replace('_pseudobulk_cn', '') for c in cn_clone.set_index(['chr', 'start', 'end']).columns.tolist()]

    # make a table of copy number breakpoints
    bk_clone = make_bk_clones(cn_clone, htert_clones)

    # compute the relative RT and CN values for each dataset
    cn_clone, rt_clone = compute_relative_rt_and_cn_clones(cn_clone, rt_clone, htert_clones)

    # merge the cn and rt tables into one long-form dataframe
    df_clone = merge_cn_and_rt_info_clones(cn_clone, rt_clone, bk_clone, htert_clones)

    # create column to denote whether a particular bin has a gain, loss, or no cna
    df_clone['cna_type'] = 'unaltered'
    for i, row in df_clone.iterrows():
        if row['relative_cn'] > 0:
            df_clone.loc[i, 'cna_type'] = 'gain'
        elif row['relative_cn'] < 0:
            df_clone.loc[i, 'cna_type'] = 'loss'

    # rename bk column
    df_clone['breakpoint'] = df_clone['pseudobulk_bk'].replace({0: 'No', 1: 'Yes'})

    # plot the RT distributions
    plot_rt_distributions(df_clone, argv)

    # plot the CN and ref RT profiles
    plot_profiles(cn_clone, rt_clone, htert_clones, argv)


if __name__ == '__main__':
    main()
