import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome import refgenome
from sklearn import preprocessing
from statannot import add_stat_annotation
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='list of pseudobulk RT profiles for each sample')
    p.add_argument('-srp', '--sample_rt_profiles', type=str, help='plot of the RT profiles for each sample')
    p.add_argument('-srd', '--sample_rt_diffs', type=str, help='violin plots of difference in RT between autosomes and chrX for samples with and without chrX CNAs')
    p.add_argument('-crp', '--clone_rt_profiles', type=str, help='plot of the RT profiles for each clone, one subpanel per sample')
    p.add_argument('-crd4', '--clone_rt_diffs_SA1054', type=str, help='violin plots of difference in RT between autosomes and chrX for clones with and without chrX CNAs for SA1054')
    p.add_argument('-crd5', '--clone_rt_diffs_SA1055', type=str, help='violin plots of difference in RT between autosomes and chrX for clones with and without chrX CNAs for SA1055')


    return p.parse_args()


color_reference = {0:'#3182BD', 1:'#9ECAE1', 2:'#CCCCCC', 3:'#FDCC8A', 4:'#FC8D59', 5:'#E34A33', 6:'#B30000', 7:'#980043', 8:'#DD1C77', 9:'#DF65B0', 10:'#C994C7', 11:'#D4B9DA'}

def get_cn_cmap(cn_data):
    min_cn = int(cn_data.min())
    max_cn = int(cn_data.max())
    assert min_cn - cn_data.min() == 0
    assert max_cn - cn_data.max() == 0
    color_list = []
    for cn in range(min_cn, max_cn+1):
        if cn > max(color_reference.keys()):
            cn = max(color_reference.keys())
        color_list.append(color_reference[cn])
    return ListedColormap(color_list)


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


def violins_with_pvals(df, x, y, hue, ax, box_pairs, order=None, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, order=order)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test, order=order,
                        text_format=text_format, loc=loc, verbose=verbose)

def plot_rt_diff_vs_chrX(df, ax, x='chrX', y='RT_diff', test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' 
    Plot the difference in pseudobulk RT values between SA039 and OV2295
    where the data is split by chrX vs autosomes
    '''
    x = x
    y = y
    hue = None
    box_pairs = [
        ('autosomes', 'chrX'),
    ]
    order = ['autosomes', 'chrX']
    violins_with_pvals(df, x, y, hue, ax, box_pairs, test=test, order=order,
                       text_format=text_format, loc=loc, verbose=verbose)
    return ax


def plot_sample_rt_profiles(rt, rt_coi, argv):
    fig, ax = plt.subplots(2, 1, figsize=(16,8), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_coi):
        dataset_id = col.split('_')[0]
        
        # plot pseudobulk rt values for this dataset
        plot_cell_cn_profile2(
            ax[0], rt, col, color='C{}'.format(i), 
            max_cn=None, scale_data=False, lines=True, label=dataset_id
        )
        # plot pseudobulk cn values for this dataset
        plot_cell_cn_profile2(
            ax[1], rt, col, color='C{}'.format(i), chromosome='X',
            max_cn=None, scale_data=False, lines=True, label=dataset_id
        )

    ax[0].set_title('Sample Pseudobulk RT')
    ax[1].set_title('Sample Pseudobulk RT')
    ax[0].set_ylabel('RT profile\n<--late | early-->')
    ax[1].set_ylabel('RT profile\n<--late | early-->')
    ax[0].legend(title='Sample ID', loc='upper left')
    ax[1].legend(title='Sample ID', loc='upper left')

    fig.savefig(argv.sample_rt_profiles, dpi=300, bbox_inches='tight')


def plot_sample_rt_diffs(rt_diff, rt_diff_coi, argv):
    fig, ax = plt.subplots(2, 2, figsize=(6,8), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_diff_coi):
        d = col.split('_')[0]
        plot_rt_diff_vs_chrX(rt_diff, ax[i], y=col)
        ax[i].set_xlabel(None)
        ax[i].set_ylabel('Difference in pseudobulk RT\n<--ref earlier | {} earlier -->'.format(d))
        ax[i].set_title('chrX RT in {}'.format(d))

    fig.savefig(argv.sample_rt_diffs, dpi=300, bbox_inches='tight')


def plot_clone_rt_profiles(rt, ax, dataset_id):
    # get the clone RT columns for this dataset
    clone_rt_cols = [c for c in rt.columns if c.startswith(dataset_id) and 'clone' in c]

    # plot each clone RT profile as a different color
    # only plot the X chromosome for now
    for i, col in enumerate(clone_rt_cols):
        clone_id = col.split('_')[2].split('clone')[1]
        plot_cell_cn_profile2(
            ax, rt, col, color='C{}'.format(i), chromosome='X',
            max_cn=None, scale_data=False, lines=True, label=clone_id
        )
    
    # format the plot
    ax.set_title('{} clone RT profiles'.format(dataset_id))
    ax.set_ylabel('RT profile\n<--late | early-->')
    ax.legend(title='Clone ID')


def plot_clone_rt_profiles_wrapper(rt, argv):
    # plot clone RT profiles for each dataset
    # find the number of datasets
    datasets = set([c.split('_')[0] for c in rt.columns if c.startswith('SA') and 'clone' not in c])
    n_datasets = len(datasets)

    # create a grid of subplots where each row is a dataset and the height is propotional to the number of datasets
    fig, ax = plt.subplots(n_datasets, 1, figsize=(16, 4*n_datasets), tight_layout=True)
    ax = ax.flatten()

    # loop through each dataset and plot the clone RT profiles
    for i, dataset_id in enumerate(datasets):
        plot_clone_rt_profiles(rt, ax[i], dataset_id)
    
    # save the figure
    fig.savefig(argv.clone_rt_profiles, dpi=300, bbox_inches='tight')


def plot_clone_rt_diffs_SA1054(rt, argv):
    ''' compute RT difference between the SA1054 clones based on allelic state '''
    # list of clones that are allelically balanced
    ref_rt_clones_SA1054 = [
        'SA1054_pseudobulk_cloneB_model_rep_state',
        'SA1054_pseudobulk_cloneC_model_rep_state',
    ]

    # list of clones that are allelically imbalanced
    rt_coi_clones_SA1054 = [
        'SA1054_pseudobulk_cloneA_model_rep_state',
        'SA1054_pseudobulk_cloneD_model_rep_state',
        'SA1054_pseudobulk_cloneE_model_rep_state',
    ]

    # TODO: use number of S-phase cells in each SA1054 clone to take a weighted average of clones D & E

    # subset the RT dataframe to only include the clones of interest or the ref clones
    rt_coi_SA1054 = ref_rt_clones_SA1054 + rt_coi_clones_SA1054 + ['chr', 'start', 'end']
    rt_diff_SA1054 = rt[rt_coi_SA1054]

    # denote which loci are part of chrX
    rt_diff_SA1054['chrX'] = rt_diff_SA1054['chr'].apply(lambda x: 'chrX' if x=='X' else 'autosomes')

    # use the ref clones to compute the reference RT
    # TODO: use the number of S-phase cells in each clone to take the weighted average of the ref clones
    rt_diff_SA1054['ref_rt'] = rt_diff_SA1054[ref_rt_clones_SA1054].mean(axis=1)

    # for all the clones in rt_coi_clones_SA1054, compute the RT difference between the clone and the ref_rt column
    rt_diff_columns_SA1054 = []
    for col in rt_coi_clones_SA1054:
        clone_id = col.split('_')[2].split('clone')[1]
        temp_diff_col = '{}_rt_diff'.format(clone_id)
        rt_diff_SA1054[temp_diff_col] = rt_diff_SA1054[col] - rt_diff_SA1054['ref_rt']
        rt_diff_columns_SA1054.append(temp_diff_col)

    # for every rt_diff_column, plot the RT difference vs chrX
    fig, ax = plt.subplots(1, 3, figsize=(9,4), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_diff_columns_SA1054):
        clone_id = col.split('_')[0]
        plot_rt_diff_vs_chrX(rt_diff_SA1054, ax[i], y=col)
        ax[i].set_xlabel(None)
        ax[i].set_ylabel('Difference in pseudobulk RT\n<--ref earlier | {} earlier -->'.format(clone_id))
        ax[i].set_title('chrX RT in SA1054 clone {}'.format(clone_id))

    fig.savefig(argv.clone_rt_diffs_SA1054, dpi=300, bbox_inches='tight')


def plot_clone_rt_diffs_SA1055(rt, argv):
    ''' compute RT difference between the SA1055 clones based on allelic state '''
    # list of clones that are allelically balanced
    ref_rt_clones_SA1055 = [
        'SA1055_pseudobulk_cloneE_model_rep_state',
        'SA1055_pseudobulk_cloneF_model_rep_state',
        'SA1055_pseudobulk_cloneG_model_rep_state',
        'SA1055_pseudobulk_cloneH_model_rep_state',
        'SA1055_pseudobulk_cloneI_model_rep_state',
    ]

    # list of clones that are allelically imbalanced
    rt_coi_clones_SA1055 = [
        'SA1055_pseudobulk_cloneA_model_rep_state',
        'SA1055_pseudobulk_cloneB_model_rep_state',
        'SA1055_pseudobulk_cloneC_model_rep_state',
        'SA1055_pseudobulk_cloneD_model_rep_state',
    ]

    # subset the RT dataframe to only include the clones of interest or the ref clones
    rt_coi_SA1055 = ref_rt_clones_SA1055 + rt_coi_clones_SA1055 + ['chr', 'start', 'end']
    rt_diff_SA1055 = rt[rt_coi_SA1055]

    # denote which loci are part of chrX
    rt_diff_SA1055['chrX'] = rt_diff_SA1055['chr'].apply(lambda x: 'chrX' if x=='X' else 'autosomes')

    # use the ref clones to compute the reference RT
    rt_diff_SA1055['ref_rt'] = rt_diff_SA1055[ref_rt_clones_SA1055].mean(axis=1)

    # for all the clones in rt_coi_clones_SA1055, compute the RT difference between the clone and the ref_rt column
    rt_diff_columns_SA1055 = []
    for col in rt_coi_clones_SA1055:
        clone_id = col.split('_')[2].split('clone')[1]
        temp_diff_col = '{}_rt_diff'.format(clone_id)
        rt_diff_SA1055[temp_diff_col] = rt_diff_SA1055[col] - rt_diff_SA1055['ref_rt']
        rt_diff_columns_SA1055.append(temp_diff_col)

    # for every rt_diff_column, plot the RT difference vs chrX
    fig, ax = plt.subplots(1, 4, figsize=(12,4), tight_layout=True)
    ax = ax.flatten()

    for i, col in enumerate(rt_diff_columns_SA1055):
        clone_id = col.split('_')[0]
        plot_rt_diff_vs_chrX(rt_diff_SA1055, ax[i], y=col)
        ax[i].set_xlabel(None)
        ax[i].set_ylabel('Difference in pseudobulk RT\n<--ref earlier | {} earlier -->'.format(clone_id))
        ax[i].set_title('chrX RT in SA1055 clone {}'.format(clone_id))

    fig.savefig(argv.clone_rt_diffs_SA1055, dpi=300, bbox_inches='tight')

def main():
    argv = get_args()

    # read in the pseudobulk RT profiles for each sample
    rt = pd.DataFrame()
    for path in argv.input:
        temp_rt = pd.read_csv(path, sep='\t')
        cols = [c for c in temp_rt.columns if c.startswith('pseudobulk') and c.endswith('model_rep_state')]
        cols.append('chr')
        cols.append('start')
        temp_rt = temp_rt[cols]
        
        # add dataset name as prefix to RT columns
        d = path.split('/')[2]
        for c in temp_rt.columns:
            if c.startswith('pseudobulk') and c.endswith('model_rep_state'):
                temp_rt.rename(columns={c: '{}_{}'.format(d, c)}, inplace=True)
        
        if rt.empty:
            rt = temp_rt
        else:
            rt = pd.merge(rt, temp_rt)

    # set chr column to category
    rt.chr = rt.chr.astype('str')
    rt.chr = rt.chr.astype('category')

    # add end position as it's necessary for plotting functions
    rt['end'] = rt['start'] + 500000 - 1

    # columns of interest for sample RT profiles
    rt_coi = [c for c in rt.columns if (c.endswith('_pseudobulk_model_rep_state'))]

    # plot the sample RT pseudobulk profiles
    plot_sample_rt_profiles(rt, rt_coi, argv)

    # samples with no CNAs on chrX should be used as a reference
    ref_rt_samples = [
        'SA039_pseudobulk_model_rep_state',
        'SA906a_pseudobulk_model_rep_state',
        'SA906b_pseudobulk_model_rep_state',
        'SA1292_pseudobulk_model_rep_state',
        'SA1188_pseudobulk_model_rep_state',
    ]

    # samples with CNAs on chrX should be used as a test
    test_rt_samples = [
        'SA1056_pseudobulk_model_rep_state',
        'SA1054_pseudobulk_model_rep_state',
        'SA1055_pseudobulk_model_rep_state',
        'OV2295_pseudobulk_model_rep_state'
    ]

    # create a new dataframe with the difference in RT values between the reference samples
    rt_diff = rt[[c for c in rt.columns if (c in rt_coi or c in ['chr', 'start', 'end'])]]

    # use the average of all the reference samples as the reference RT profile
    rt_diff['ref_rt'] = rt[ref_rt_samples].mean(axis=1)

    # denote which loci are on chrX
    rt_diff['chrX'] = rt_diff['chr'].apply(lambda x: 'chrX' if x=='X' else 'autosomes')

    # calculate the difference in RT values between the reference and test samples
    for col in test_rt_samples:
        dataset_id = col.split('_')[0]
        rt_diff['{}_rt_diff'.format(dataset_id)] = rt_diff[col] - rt_diff['ref_rt']
    
    # denote which columns are of interest for plotting
    rt_diff_coi = [
        'SA1056_rt_diff',
        'SA1054_rt_diff',
        'SA1055_rt_diff',
        'OV2295_rt_diff'
    ]

    # plot the difference in RT values between the reference and test samples
    plot_sample_rt_diffs(rt_diff, rt_diff_coi, argv)

    # plot the clone RT profiles for chrX
    plot_clone_rt_profiles_wrapper(rt, argv)

    # plot clone-specific RT differences on chrX for SA1054 and SA1055
    plot_clone_rt_diffs_SA1054(rt, argv)
    plot_clone_rt_diffs_SA1055(rt, argv)


if __name__ == '__main__':
    main()
