import pandas as pd
import numpy as np
from scipy.stats import linregress
from argparse import ArgumentParser
import matplotlib.pyplot as plt
# save svg text as text, not curves
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams["font.family"] = 'Arial'
plt.rcParams['pdf.use14corefonts'] = True
import seaborn as sns


def get_args():
    p = ArgumentParser()

    p.add_argument('-sb', '--sample_bafs', type=str, help='Mean B-allele frequencies for DNA and RNA per sample and chromsome arm')
    p.add_argument('-cb', '--clone_bafs', type=str, help='Mean B-allele frequencies for DNA per clone and chromsome arm')
    p.add_argument('-sr', '--sample_rt', type=str, help='chrX RT delays for sample profiles')
    p.add_argument('-cr', '--clone_rt', type=str, help='chrX RT delays for clone profiles')
    p.add_argument('-c', '--counts', type=str, help='merged counts for all samples')
    p.add_argument('-o', '--output', type=str, help='output figure path')

    return p.parse_args()


def merge_sample_rt_and_bafs(sample_mean_rt, sample_bafs):
    ''' Query and pivot the dataframes such that each row has both the chrX RT delay and BAFs for a given sample '''
    # subset df to just chr=='X'
    sample_df = sample_bafs.query('chr=="X"')[['patient', 'chr_arm', 'dna_baf_mean', 'rna_baf_mean', 'dna_baf_std', 'rna_baf_std']]
    # pivot such that chrXp_dna_baf, chrXq_dna_baf, chrXp_rna_baf, chrXq_rna_baf are columns
    sample_df = sample_df.pivot(index='patient', columns='chr_arm', values=['dna_baf_mean', 'rna_baf_mean', 'dna_baf_std', 'rna_baf_std'])
    # flatten the column names
    sample_df.columns = ['_'.join(col).strip() for col in sample_df.columns.values]
    # reset the index
    sample_df = sample_df.reset_index().rename(columns={'patient': 'dataset'})
    # merge sample_df with sample_mean_rt
    sample_df = pd.merge(sample_df, sample_mean_rt)

    return sample_df


def merge_sample_mapping(sample_df):
    # create a mapping of sample_id to cell/tumor type, condition, and signature
    sample_id_to_type = {
        'SA1047': 'HGSOC',
        'SA1049': 'HGSOC',
        'SA1050': 'HGSOC',
        'SA1051': 'HGSOC',
        'SA1052': 'HGSOC',
        'SA1053': 'HGSOC',
        'SA1091': 'HGSOC',
        'SA1093': 'HGSOC',
        'SA1096': 'HGSOC',
        'SA1162': 'HGSOC',
        'SA1181': 'HGSOC',
        'SA1182': 'HGSOC',
        'SA1184': 'HGSOC',
        'SA501': 'TNBC',
        'SA530': 'TNBC',
        'SA604': 'TNBC',
        'SA039': 'hTERT',
        'SA906a': 'hTERT',
        'SA906b': 'hTERT',
        'SA1188': 'hTERT',
        'SA1292': 'hTERT',
        'SA1054': 'hTERT',
        'SA1055': 'hTERT',
        'SA1056': 'hTERT',
        'OV2295': 'OV2295',
        'SA535': 'TNBC',
        'SA609': 'TNBC',
        'SA1035': 'TNBC',
        'T47D': 'T47D',
        'GM18507': 'GM18507'
    }
    sample_id_to_condition = {
        'SA1047': 'PDX',
        'SA1049': 'PDX',
        'SA1050': 'PDX',
        'SA1051': 'PDX',
        'SA1052': 'PDX',
        'SA1053': 'PDX',
        'SA1091': 'PDX',
        'SA1093': 'PDX',
        'SA1096': 'PDX',
        'SA1162': 'PDX',
        'SA1181': 'PDX',
        'SA1182': 'PDX',
        'SA1184': 'PDX',
        'SA501': 'PDX',
        'SA530': 'PDX',
        'SA604': 'PDX',
        'SA039': 'Line',
        'SA906a': 'Line',
        'SA906b': 'Line',
        'SA1188': 'Line',
        'SA1292': 'Line',
        'SA1054': 'Line',
        'SA1055': 'Line',
        'SA1056': 'Line',
        'OV2295': 'Line',
        'SA535': 'PDX',
        'SA609': 'PDX',
        'SA1035': 'PDX',
        'T47D': 'Line',
        'GM18507': 'Line'
    }
    sample_id_to_signature = {
        'SA1047': 'TD',
        'SA1049': 'FBI',
        'SA1050': 'HRD-Dup',
        'SA1051': 'HRD-Dup',
        'SA1052': 'HRD-Dup',
        'SA1053': 'HRD-Dup',
        'SA1091': 'FBI',
        'SA1093': 'TD',
        'SA1096': 'FBI',
        'SA1162': 'FBI',
        'SA1181': 'HRD-Dup',
        'SA1182': 'FBI',
        'SA1184': 'HRD-Dup',
        'SA501': 'HRD-Dup',
        'SA530': 'FBI',
        'SA604': 'FBI',
        'SA039': 'N/A',
        'SA906a': 'N/A',
        'SA906b': 'N/A',
        'SA1188': 'N/A',
        'SA1292': 'N/A',
        'SA1054': 'N/A',
        'SA1055': 'N/A',
        'SA1056': 'N/A',
        'OV2295': 'N/A',
        'SA535': 'HRD-Dup',
        'SA609': 'FBI',
        'SA1035': 'N/A',
        'T47D': 'N/A',
        'GM18507': 'N/A'
    }

    sample_mapping = pd.DataFrame(sample_id_to_condition.items(), columns=['dataset', 'condition'])
    sample_mapping['type'] = sample_mapping['dataset'].map(sample_id_to_type)
    sample_mapping['signature'] = sample_mapping['dataset'].map(sample_id_to_signature)
    
    # merge sample_df with sample_mapping
    sample_df = pd.merge(sample_df, sample_mapping)

    return sample_df, sample_mapping


def plot_hTERT_sample_bafs_vs_rt(sample_df, ax):
    # plot the hTERT samples
    sns.regplot(y='mean_chrX_rt_delay', x='dna_chrX_baf_mean', data=sample_df.query('type=="hTERT"'), ax=ax, ci=None)
    ax.set_ylabel('chrX RT delay\n<--delayed | advanced -->')
    ax.set_xlabel('chrX DNA BAF')
    ax.set_title('hTERT samples')
    # add r and p-value annotations in the bottom-left corner
    reg = linregress(sample_df.query('type=="hTERT"')['dna_chrX_baf_mean'].values, sample_df.query('type=="hTERT"')['mean_chrX_rt_delay'].values)
    ax.text(0.05, 0.05, 'r={:.2f}\np={:.2e}'.format(reg.rvalue, reg.pvalue), transform=ax.transAxes)
    # edit the x-and y-axis limits to give 0.05 padding on all sides compared to the current limits
    ax.set_xlim([ax.get_xlim()[0] - 0.05, ax.get_xlim()[1] + 0.05])
    ax.set_ylim([max(-.4, ax.get_ylim()[0] - 0.05), ax.get_ylim()[1] + 0.05])


def plot_SA1054_subclonal_bafs_vs_rt(clone_df, ax):
    # scatterplot of the data
    sns.scatterplot(y='mean_chrX_rt_delay', x='dna_chrX_baf_mean', hue='clone_id', data=clone_df.query('dataset=="SA1054"'), ax=ax, size='num_cells_s', sizes=(10, 100))

    sns.regplot(y='mean_chrX_rt_delay', x='dna_chrX_baf_mean', data=clone_df.query('dataset=="SA1054"'), scatter=False, ax=ax)

    # fit a linear regression model to the data
    res = linregress(clone_df.query('dataset=="SA1054"')['dna_chrX_baf_mean'].values, clone_df.query('dataset=="SA1054"')['mean_chrX_rt_delay'].values)
    # add r- and p-value annotations in the bottom-left corner
    ax.text(0.05, 0.05, 'r={:.2f}\np={:.2e}'.format(res.rvalue, res.pvalue), transform=ax.transAxes)

    # edit the x-and y-axis limits to give 0.05 padding on all sides compared to the current limits
    ax.set_xlim([ax.get_xlim()[0] - 0.05, ax.get_xlim()[1] + 0.05])
    ax.set_ylim([ax.get_ylim()[0] - 0.05, ax.get_ylim()[1] + 0.05])

    ax.set_ylabel('chrX RT delay\n<--delayed | advanced-->')
    ax.set_xlabel('chrX DNA BAF')
    ax.set_title('SA1054 clones')

    # add the legend to the right of the plot
    ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)


def plot_sample_rt_delays(sample_df, ax):
    # make a barplot of 'mean_chrX_rt_delay' on the y-axis and 'dataset' on the x-axis with the bars colored by 'type'
    sns.barplot(y='mean_chrX_rt_delay', x='dataset', hue='type', data=sample_df.sort_values(by=['type', 'dataset']), ax=ax)
    ax.set_ylabel('chrX RT delay\n<--delayed | advanced -->')
    ax.set_xlabel('Sample ID')
    ax.set_title('All samples')
    # rotate the x-tick labels by 45 degrees
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')


def main():
    argv = get_args()

    # load BAFs, RT delays, and counts
    sample_bafs = pd.read_csv(argv.sample_bafs)
    clone_bafs = pd.read_csv(argv.clone_bafs)
    sample_mean_rt = pd.read_csv(argv.sample_rt)
    clone_mean_rt = pd.read_csv(argv.clone_rt)
    counts = pd.read_csv(argv.counts)

    # merge the clone RT delays and BAFs
    clone_df = clone_bafs.query('chr=="X"')[['dataset', 'clone_id', 'dna_baf_mean']].rename(columns={
        'dna_baf_mean': 'dna_chrX_baf_mean'
    }).reset_index(drop=True)
    clone_df = pd.merge(clone_df, clone_mean_rt)

    # merge sample RT delays and BAFs
    sample_df = merge_sample_rt_and_bafs(sample_mean_rt, sample_bafs)
    # merge signature, cell type and condition info with sample ID
    sample_df, sample_mapping = merge_sample_mapping(sample_df)
    # merge the same sample info and counts into clone_df
    clone_df = pd.merge(clone_df, sample_mapping)
    clone_df = pd.merge(clone_df, counts, how='left')

    # create an 8.5 x 11 figure with 4 rows and 2 columns
    fig, axes = plt.subplots(4, 2, figsize=(8.5, 11), tight_layout=True)
    # merge the two subplots in the 2nd row into one subplot
    axes[1, 0].remove()
    axes[1, 1].remove()
    axes[1, 0] = fig.add_subplot(4, 2, (3, 4))

    # the top left plot is the chrX RT delay vs DNA BAF for all hTERT sample pseudobulks
    plot_hTERT_sample_bafs_vs_rt(sample_df, axes[0, 0])

    # the top right plot is the SA1054 clone RT delay vs DNA BAF
    # plot_SA1054_subclonal_bafs_vs_rt(clone_df, axes[0, 1])

    # create a barplot of all the sample RT delays, sorted and colored by type
    plot_sample_rt_delays(sample_df, axes[1, 0])

    # TODO: add the remaining subplots


    # save the figure
    fig.savefig(argv.output, dpi=300, bbox_inches='tight')




if __name__ == '__main__':
    main()