import matplotlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.plot_utils import plot_cell_cn_profile2
from common.colors import get_cna_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('-ir', '--input_rt', type=str, nargs='+', help='pseduobulk rt profiles for each dataset')
    p.add_argument('-ic', '--input_cn', type=str, nargs='+', help='pseduobulk cn profiles for each dataset')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='list of dataset names')
    p.add_argument('--plot1', help='plots showing the RT distribution against CNA types and breakpoints')
    p.add_argument('--plot2', help='plots reference RT and pseudobulk CN profiles')

    return p.parse_args()


def load_rt_data(argv):
    # load dataset pseudobulk rt profiles
    rt = pd.DataFrame()
    for rt_path, d in zip(argv.input_rt, argv.datasets):
        temp_rt = pd.read_csv(rt_path, sep='\t')
        if 'SA039' in rt_path:
            temp_rt = temp_rt[[
                'chr', 'start', 'pseudobulk_cloneA_model_rep_state', 'pseudobulk_cloneA_hours', 'pseudobulk_model_rep_state', 'pseudobulk_hours'
            ]]
            temp_rt.rename(columns={
                'pseudobulk_cloneA_model_rep_state': '{}_cloneA_pseudobulk_rt'.format(d),
                'pseudobulk_cloneA_hours': '{}_cloneA_pseudobulk_rt_hours'.format(d)
            }, inplace=True)
        else:
            temp_rt = temp_rt[[
                'chr', 'start', 'pseudobulk_model_rep_state', 'pseudobulk_hours',
            ]]
        temp_rt.rename(columns={
            'pseudobulk_model_rep_state': '{}_pseudobulk_rt'.format(d),
            'pseudobulk_hours': '{}_pseudobulk_rt_hours'.format(d)
        }, inplace=True)
        
        if rt.empty:
            rt = temp_rt
        else:
            rt = pd.merge(rt, temp_rt)

    # set chr column to category
    rt.chr = rt.chr.astype('str')
    rt.chr = rt.chr.astype('category')

    rt['end'] = rt['start'] + 500000 - 1
            
    return rt
    

def load_cn_data(argv):
    # load dataset pseudobulk cn profiles
    cn = pd.DataFrame()

    col_name_change = {}
    for cn_path, d in zip(argv.input_cn, argv.datasets):
        temp_cn = pd.read_csv(cn_path, sep='\t')
        temp_cn = temp_cn.set_index(['chr', 'start'])
        
        # find the sample level pseudobulk cn column
        # note that it can start with dataset or sample prefix
        col_of_interest = [c for c in temp_cn.columns if c.startswith('dataset')]
        if len(col_of_interest)==0:
            col_of_interest = [c for c in temp_cn.columns if c.startswith('sample')]
        temp_cn = temp_cn[col_of_interest].reset_index()
        
        col_name_change[col_of_interest[0]] = '{}_pseudobulk_cn'.format(d)
        
        if cn.empty:
            cn = temp_cn
        else:
            cn = pd.merge(cn, temp_cn)


    # set chr column to category
    cn.chr = cn.chr.astype('str')
    cn.chr = cn.chr.astype('category')

    cn['end'] = cn['start'] + 500000 - 1

    cn = cn.rename(columns=col_name_change)

    return cn


def make_bk(cn, argv):
    # create a table of copy number breakpoints for sample pseudobulk data
    # copy over the chromosome, start, and end columns
    cn_breakpoints = cn[['chr', 'start', 'end']].copy()
    for d in argv.datasets:
        # create a column for the copy number breakpoints
        cn_breakpoints['{}_pseudobulk_cn_breakpoints'.format(d)] = cn['{}_pseudobulk_cn'.format(d)].diff().fillna(0).astype('int')
        # convert all nonzero values to 1
        cn_breakpoints['{}_pseudobulk_cn_breakpoints'.format(d)] = cn_breakpoints['{}_pseudobulk_cn_breakpoints'.format(d)].apply(lambda x: 1 if x!=0 else 0)

    return cn_breakpoints


def compute_relative_rt_and_cn(cn, rt, argv):
    for d in argv.datasets:
        ref_rt_col = 'SA039_cloneA_pseudobulk_rt'
        ref_cn_col = 'SA039_pseudobulk_cn'
        
        temp_rt_col = '{}_pseudobulk_rt'.format(d)
        temp_cn_col = '{}_pseudobulk_cn'.format(d)
        
        relative_rt_col = '{}_pseudobulk_relative_rt'.format(d)
        relative_cn_col = '{}_pseudobulk_relative_cn'.format(d)
        
        rt[relative_rt_col] = rt[temp_rt_col] - rt[ref_rt_col]
        cn[relative_cn_col] = (cn[temp_cn_col] / cn[temp_cn_col].mode().values[0]) - (cn[ref_cn_col] / cn[ref_cn_col].mode().values[0])
    
    return cn, rt


def merge_cn_and_rt_info(cn, rt, bk, argv):
    # merge the cn and rt tables into one long-form dataframe
    df = []

    for i, d in enumerate(argv.datasets):
        ref_rt_col = 'SA039_cloneA_pseudobulk_rt'
        ref_cn_col = 'SA039_pseudobulk_cn'
        
        temp_rt_col = '{}_pseudobulk_rt'.format(d)
        temp_cn_col = '{}_pseudobulk_cn'.format(d)
        temp_bk_col = '{}_pseudobulk_cn_breakpoints'.format(d)
        
        relative_rt_col = '{}_pseudobulk_relative_rt'.format(d)
        relative_cn_col = '{}_pseudobulk_relative_cn'.format(d)
        
        temp_df = pd.DataFrame({
            'chr': cn['chr'], 'start': cn['start'], 'end': cn['end'], 'dataset': [d]*cn.shape[0],
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
        ('loss', 'neutral'),
        ('neutral', 'gain'),
    ]
    order = ['loss', 'neutral', 'gain']
    violins_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose, palette=cna_cmap)
    return ax


def plot_WT_rt_vs_bk(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' 
    Plot the distribution of SA039 (WT) RT values against the ensuing breakpoints that emerge
    in the other hTERT cell lines. This data will show whether CN breakpoints preferentially
    emerge in early or late replicating regions.
    '''
    x = "breakpoint"
    y = "WT_pseudobulk_rt"
    hue = None
    box_pairs = [
        ('No', 'Yes'),
    ]
    order = ['No', 'Yes']
    violins_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose)
    return ax



def plot_rt_distributions(df, argv):
    fig, ax = plt.subplots(2, 2, figsize=(8,8), tight_layout=True)
    ax = ax.flatten()

    # violin plots of RT distributions split by CNA type
    plot_WT_rt_vs_cna_type(df.query('dataset!="SA039"'), ax[0])
    ax[0].set_xlabel('Clonal CN in hTERT mutant cell lines')
    ax[0].set_ylabel('Ancestral hTERT RT\n<--late | early-->')
    ax[0].set_title('Location of clonal CNAs')

    # histogram of RT values split by CNA type
    df['CNA type'] = df['cna_type']
    sns.histplot(
        data=df.query('dataset!="SA039"'), x='WT_pseudobulk_rt', hue='CNA type', 
        common_norm=False, stat='density', ax=ax[1],
        palette=get_cna_cmap()
    )
    ax[1].set_xlabel('Ancestral hTERT RT\n<--late | early-->')
    ax[1].set_title('Location of clonal CNAs')

    # violin plots of RT distributions split by breakpoint presence
    plot_WT_rt_vs_bk(df.query('dataset!="SA039"'), ax[2])
    ax[2].set_xlabel('Clonal CN breakpoint in hTERT mutant cell lines')
    ax[2].set_ylabel('Ancestral hTERT RT\n<--late | early-->')
    ax[2].set_title('Location of clonal CNA breakpoints')

    # histogram of RT values split by breakpoint presence
    sns.histplot(
        data=df.query('dataset!="SA039"'), x='WT_pseudobulk_rt', hue='breakpoint', 
        common_norm=False, stat='density', ax=ax[3],
        palette={'No': 'C0', 'Yes': 'C1'}
    )
    ax[3].set_xlabel('Ancestral hTERT RT\n<--late | early-->')
    ax[3].set_title('Location of clonal CNA breakpoints')

    # save the figure
    fig.savefig(argv.plot1, dpi=300, bbox_inches='tight')


def plot_profiles(cn, rt, argv):
    # plot the reference RT profile and the relative CN profiles for each dataset
    # get the reference pseudobulk rt column
    ref_rt_col = 'SA039_cloneA_pseudobulk_rt'

    fig, ax = plt.subplots(2, 1, figsize=(16,8), tight_layout=True)
    ax = ax.flatten()

    # plot pseudobulk rt values of the reference WT dataset
    plot_cell_cn_profile2(
        ax[0], rt, ref_rt_col,
        max_cn=None, scale_data=False, lines=True
    )

    for i, d in enumerate(argv.datasets):
        relative_cn_col = '{}_pseudobulk_relative_cn'.format(d)
            
        # plot pseudobulk cn values for this dataset
        plot_cell_cn_profile2(
            ax[1], cn, relative_cn_col, color='C{}'.format(i),
            max_cn=None, scale_data=False, lines=True, label=d
        )

    ax[0].set_title('Ancestral hTERT RT (SA039 clone A)')
    ax[1].set_title('Sample Pseudobulk CN')
    ax[0].set_ylabel('Pseudobulk RT\n<--late | early-->')
    ax[1].set_ylabel('CN relative to hTERT WT and ploidy\n<--loss | gain-->')
    ax[1].legend(title='Sample ID')

    # manually set the y-ticks for ax[1] to range from -1 to 1, spaced by 0.2
    ax[1].set_yticks(np.arange(-1, 1.2, 0.2))
    ax[1].spines['left'].set_bounds(-1, 1)

    fig.savefig(argv.plot2, dpi=300, bbox_inches='tight')


def main():
    argv = get_args()

    # load the rt and cn data
    rt = load_rt_data(argv)
    cn = load_cn_data(argv)

    # make a table of copy number breakpoints
    bk = make_bk(cn, argv)

    # compute the relative RT and CN values for each dataset
    cn, rt = compute_relative_rt_and_cn(cn, rt, argv)

    # merge the cn and rt tables into one long-form dataframe
    df = merge_cn_and_rt_info(cn, rt, bk, argv)

    # create column to denote whether a particular bin has a gain, loss, or no cna
    df['cna_type'] = 'neutral'
    for i, row in df.iterrows():
        if row['relative_cn'] > 0:
            df.loc[i, 'cna_type'] = 'gain'
        elif row['relative_cn'] < 0:
            df.loc[i, 'cna_type'] = 'loss'

    # rename bk column
    df['breakpoint'] = df['pseudobulk_bk'].replace({0: 'No', 1: 'Yes'})

    # plot the RT distributions
    plot_rt_distributions(df, argv)

    # plot the CN and ref RT profiles
    plot_profiles(cn, rt, argv)


if __name__ == '__main__':
    main()
