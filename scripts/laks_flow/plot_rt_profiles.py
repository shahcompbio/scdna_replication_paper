from argparse import ArgumentParser
import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_cell_cn_profile
from matplotlib.colors import ListedColormap
from scgenome import refgenome
from sklearn import preprocessing


def get_args():
    p = ArgumentParser()

    p.add_argument('rt_bulks', help='input rt pseudobulks from multiple cell lines and model versions')
    p.add_argument('rt_diff_split', help='difference in rt values between the two cell lines for the split pyro model')
    p.add_argument('rt_diff_joint', help='difference in rt values between the two cell lines for the joint pyro model')
    p.add_argument('rt_corr', help='pairwise correlation between all rt profiles')
    p.add_argument('rt_split_chr1', help='both cell line chr1 RT profiles for split model')
    p.add_argument('rt_joint_chr1', help='both cell line chr1 RT profiles for joint model')


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


def plot_chr1_by_cell_line(df, argv):
    # plot chr1 bulk RT for each cell line inferred with the joint model
    fig, ax = plt.subplots(1, 1, figsize=(16,4))
    plot_cell_cn_profile2(ax, df, 'rt_joint_T47D', color='C0', label='joint T47D', max_cn=None, scale_data=False, lines=True, chromosome='1')
    plot_cell_cn_profile2(ax, df, 'rt_joint_GM18507', color='C1', label='joint GM18507', max_cn=None, scale_data=False, lines=True, chromosome='1')
    ax.set_ylabel('Pseudobulk RT')
    ax.set_title('Joint Model')
    ax.legend()
    fig.savefig(argv.rt_joint_chr1, bbox_inches='tight', dpi=300)

    # plot chr1 bulk RT for each cell line inferred with the split model
    fig, ax = plt.subplots(1, 1, figsize=(16,4))
    plot_cell_cn_profile2(ax, df, 'rt_split_T47D', color='C0', label='split T47D', max_cn=None, scale_data=False, lines=True, chromosome='1')
    plot_cell_cn_profile2(ax, df, 'rt_split_GM18507', color='C1', label='split GM18507', max_cn=None, scale_data=False, lines=True, chromosome='1')
    ax.set_ylabel('Pseudobulk RT')
    ax.set_title('Split Model')
    ax.legend()
    fig.savefig(argv.rt_split_chr1, bbox_inches='tight', dpi=300)


def plot_rt_corr(df, argv):
    """ plot a heatmap of all the pairwise RT correlations for each cell line when the joint and split methods are used """
    # subset the dataframe to only the columns that should be plotted
    df = df[['gc', 'rt_joint_T47D', 'rt_split_T47D', 'rt_joint_GM18507', 'rt_split_GM18507']]
    # rename the columns to be more readable
    df.columns = ['gc', 'T47D joint', 'T47D split', 'GM18507 joint', 'GM18507 split']
    # calculate the Pearson correlation matrix
    corr = df.corr()

    # mask the upper triangle of the heatmap
    mask = np.zeros_like(corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    fig = plt.figure(figsize=(6, 6))

    # plot the heatmap, including the first 2 digits of the correlation values, the mask, and whitespace between the cells
    sns.heatmap(corr, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='.2f')
    plt.title('Pseudobulk RT correlations\nFlow sorted data passed into PERT joint vs split by cell line')

    # rotate the y-tick labels to read from left to right
    plt.yticks(rotation=0)

    fig.savefig(argv.rt_corr, bbox_inches='tight', dpi=300)


def plot_rt_diff(df, argv):
    # plot the difference in RT values between the two cell lines
    # for the joint model
    fig, ax = plt.subplots(1, 1, figsize=(16,4))
    plot_cell_cn_profile2(ax, df, 'rt_diff_joint', color='C3', max_cn=None, scale_data=False, lines=True)
    ax.set_ylabel('RT Diff\n<--GM18507 earlier | T47D earlier -->')
    ax.set_title('Joint Model')
    fig.savefig(argv.rt_diff_joint, bbox_inches='tight', dpi=300)

    # for the split model
    fig, ax = plt.subplots(1, 1, figsize=(16,4))
    plot_cell_cn_profile2(ax, df, 'rt_diff_split', color='C3', max_cn=None, scale_data=False, lines=True)
    ax.set_ylabel('RT Diff\n<--GM18507 earlier | T47D earlier -->')
    ax.set_title('Split Model')
    fig.savefig(argv.rt_diff_split, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()
    df = pd.read_csv(argv.rt_bulks, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # create separate plots
    plot_rt_corr(df, argv)
    plot_rt_diff(df, argv)
    plot_chr1_by_cell_line(df, argv)


if __name__ == '__main__':
    main()
