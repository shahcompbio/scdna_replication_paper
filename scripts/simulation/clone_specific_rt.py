import matplotlib
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome import refgenome
from sklearn import preprocessing
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('P10_input', type=str, help='pseudobulk RT profiles for simulated dataset P10')
    p.add_argument('P11_input', type=str, help='pseudobulk RT profiles for simulated dataset P11')
    p.add_argument('output_corr', type=str, help='Heatmap of pairwise correlations between true and inferred clone RT profiles')
    p.add_argument('output_profiles', type=str, help='Plot the clone RT profiles directly')

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


def plot_clone_rt_profiles(df, argv):
    # plot clone RT profiles for each dataset
    # find the number of datasets
    inferred_cols = [c for c in df.columns if 'inferred' in c]
    true_cols = [c for c in df.columns if 'true' in c]

    # create a grid of subplots where each row is a dataset and the height is propotional to the number of datasets
    fig, ax = plt.subplots(2, 2, figsize=(32, 8), tight_layout=True)
    ax = ax.flatten()

    # plot the inferred RT profiles for each dataset
    for i, col in enumerate(inferred_cols):
        label = col.replace('_', ' ')
        # plot the inferred RT profile for the whole genome
        plot_cell_cn_profile2(
            ax[0], df, col, color='C{}'.format(i),
            max_cn=None, scale_data=False, lines=True, label=label
        )
        # plot the inferred RT profile for the chr1
        plot_cell_cn_profile2(
            ax[1], df, col, color='C{}'.format(i), chromosome='1',
            max_cn=None, scale_data=False, lines=True, label=label
        )
    
    # plot the true RT profiles for each dataset
    for i, col in enumerate(true_cols):
        label = col.replace('_', ' ')
        # plot the inferred RT profile for the whole genome
        plot_cell_cn_profile2(
            ax[2], df, col, color='C{}'.format(i),
            max_cn=None, scale_data=False, lines=True, label=label
        )
        # plot the inferred RT profile for the chr1
        plot_cell_cn_profile2(
            ax[3], df, col, color='C{}'.format(i), chromosome='1',
            max_cn=None, scale_data=False, lines=True, label=label
        )
    
    for a in [ax[0], ax[1]]:
        a.set_title('Inferred RT profiles')
        a.set_ylabel('RT profile\n<--late | early-->')
        a.legend(title='Clone ID', loc='upper left')
    for a in [ax[2], ax[3]]:
        a.set_title('True RT profiles')
        a.set_ylabel('RT profile\n<--late | early-->')
        a.legend(title='Clone ID', loc='upper left')
    
    # save the figure
    fig.savefig(argv.output_profiles, dpi=300, bbox_inches='tight')


def main():
    argv = get_args()

    # read in the pseudobulk RT profiles for each dataset
    df10 = pd.read_csv(argv.P10_input, sep='\t')
    df11 = pd.read_csv(argv.P11_input, sep='\t')

    # subset to just the acceptable columns in each dataset
    good_cols = [
        'chr', 'start', 
        'pseduobulk_cloneA_model_rep_state', 'true_pseduobulk_cloneA_true_rep',
        'pseduobulk_cloneB_model_rep_state', 'true_pseduobulk_cloneB_true_rep',
        'pseduobulk_cloneC_model_rep_state', 'true_pseduobulk_cloneC_true_rep',
        'pseduobulk_cloneD_model_rep_state', 'true_pseduobulk_cloneD_true_rep'
    ]

    df10 = df10[good_cols]
    df11 = df11[good_cols]

    # remove the 'true_' from the 'true_pseduobulk_' prefix from the true replication columns
    df10.columns = df10.columns.str.replace('true_pseduobulk', 'pseduobulk')
    df11.columns = df11.columns.str.replace('true_pseduobulk', 'pseduobulk')

    # correct the pseudobulk typo
    df10.columns = df10.columns.str.replace('pseduobulk', 'pseudobulk')
    df11.columns = df11.columns.str.replace('pseduobulk', 'pseudobulk')

    # create a dict that maps the clone names to the cell lines used in the simulation
    # for dataset 10
    clone_to_cell_line_10 = {
        'cloneA': 'gm06990_P10',
        'cloneB': 'gm12801_P10',
        'cloneC': 'gm12812_P10',
        'cloneD': 'gm12813_P10'
    }

    # for dataset 11
    clone_to_cell_line_11 = {
        'cloneA': 'bj_P11',
        'cloneB': 'mcf7_P11',
        'cloneC': 'hepg2_P11',
        'cloneD': 'gm12813_P11'
    }

    # use these two dicts to rename the columns in each dataset wherever the substring 'clone' appears
    df10 = df10.rename(columns={col: col.replace('clone', clone_to_cell_line_10[col.split('_')[1]]) for col in df10.columns if 'clone' in col})
    df11 = df11.rename(columns={col: col.replace('clone', clone_to_cell_line_11[col.split('_')[1]]) for col in df11.columns if 'clone' in col})

    # merge the two datasets on the chr and start columns
    df = pd.merge(df10, df11, on=['chr', 'start'])

    # convert chr column to category
    df['chr'] = df['chr'].astype('category')
    df['end'] = df['start'] + 500000 - 1

    # drop the 'pseudobulk_' prefix from the columns
    df.columns = df.columns.str.replace('pseudobulk_', '')

    # rename the 'model_rep_state' suffix to 'inferred'
    df.columns = df.columns.str.replace('model_rep_state', 'inferred')

    # rename the 'true_rep' suffix to 'true'
    df.columns = df.columns.str.replace('true_rep', 'true')

    # plot the clone RT profiles directly
    plot_clone_rt_profiles(df, argv)

    # compute the pairwise correlation between all of the columns except the chr and start columns
    corr = df.corr()

    # drom the 'start' row and column from the correlation matrix
    corr = corr.drop('start', axis=0)
    corr = corr.drop('start', axis=1)

    # replace the underscores in the column and row names with spaces
    corr.columns = corr.columns.str.replace('_', ' ')
    corr.index = corr.index.str.replace('_', ' ')

    # plot the correlation matrix as a heatmap
    fig = plt.figure(figsize=(10, 10))

    # mask the upper triangle of the heatmap
    mask = np.zeros_like(corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    # plot the heatmap, including the first 2 digits of the correlation values
    sns.heatmap(corr, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='.2f')
    plt.title('Pseudobulk RT correlations\nSimulated datasets P10 and P11 with clone-specific RT')

    # save the figure
    fig.savefig(argv.output_corr, bbox_inches='tight', dpi=300)

if __name__ == '__main__':
    main()
