import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
from scipy.spatial.distance import squareform, pdist


def get_args():
    p = ArgumentParser()

    p.add_argument('-ic', '--input_cn', type=str, nargs='+', help='list of pseudobulk CN profiles for each sample')
    p.add_argument('-ir', '--input_rt', type=str, nargs='+', help='list of pseudobulk RT profiles for each sample')
    p.add_argument('-c', '--counts', type=str, help='cell cycle clone counts across all samples')
    p.add_argument('-sc', '--sample_corrs', type=str, help='plot of CN and RT correlations between samples')
    p.add_argument('-cc', '--clone_corrs', type=str, help='plot of CN and RT correlations between clones')
    p.add_argument('-rc', '--rx_corrs', type=str, help='plot of RT correlations between samples and treatment status')

    return p.parse_args()


def read_rt_data(argv):
    # read in the pseudobulk RT profiles for each sample
    rt = pd.DataFrame()
    for path in argv.input_rt:
        temp_rt = pd.read_csv(path, sep='\t')
        cols = [c for c in temp_rt.columns if c.startswith('pseudobulk') and ('hours' not in c)]
        cols.append('chr')
        cols.append('start')
        temp_rt = temp_rt[cols]
        
        # add dataset name as prefix to RT columns
        d = path.split('/')[2]
        for c in temp_rt.columns:
            if c.startswith('pseudobulk') and ('hours' not in c):
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

    return rt


def read_cn_data(argv):
    # read in the pseudobulk CN profiles for each sample
    cn = pd.DataFrame()

    for file in argv.input_cn:
        temp_cn = pd.read_csv(file, sep='\t')

        # find the dataset name given the path
        d = file.split('/')[2]

        # rename the 'sample_{dataset}' column to '{dataset}_sample_cn'
        if 'sample_{}'.format(d) in temp_cn.columns:
            temp_cn.rename(columns={'sample_{}'.format(d): '{}_sample_cn'.format(d)}, inplace=True)
        # sometimes the sample column is named 'dataset_{dataset}'
        elif 'dataset_{}'.format(d) in temp_cn.columns:
            temp_cn.rename(columns={'dataset_{}'.format(d): '{}_sample_cn'.format(d)}, inplace=True)

        # add the dataset name as prefix to all columns with 'clone' in the name
        cols = [c for c in temp_cn.columns if c.startswith('clone')]
        for c in cols:
            temp_cn.rename(columns={c: '{}_{}_cn'.format(d, c)}, inplace=True)
        
        # merge the temp_cn dataframe with the cn dataframe if cn is not empty
        if cn.empty:
            cn = temp_cn
        else:
            cn = pd.merge(cn, temp_cn)

    # set chr column to category
    cn.chr = cn.chr.astype('str')
    cn.chr = cn.chr.astype('category')

    cn['end'] = cn['start'] + 500000 - 1

    return cn


def plot_sample_corrs(cn, rt, sample_cn_cols, sample_rt_cols, argv):
    ''' Create a figure that shows the correlation between the sample CN and RT profiles. '''
    fig, ax = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)

    # pairwise distance between copy number profiles
    cn_dist = pd.DataFrame(
        squareform(pdist(cn.T.loc[sample_cn_cols])),
        columns = sample_cn_cols,
        index = sample_cn_cols
    )
    # subtract the max value from the distance matrix to make the heatmap more intuitive
    # this way the larger the value, the more similar the profiles
    cn_dist = cn_dist.max().max() - cn_dist
    # mask the upper triangle of the heatmap
    mask = np.zeros_like(cn_dist, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    # plot the heatmap
    sns.heatmap(cn_dist, square=False, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='1.0f', ax=ax[0])
    ax[0].set_title('Similarity in sample CN')

    # correlation between RT profiles
    rt_corrs = rt[sample_rt_cols].corr()
    mask = np.zeros_like(rt_corrs, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(rt_corrs, square=False, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='.2f', ax=ax[1])
    ax[1].set_title('Correlation in sample RT')

    # only include the sample_id prefix in the xticklabels and yticklabels
    for a in ax:
        a.set_xticklabels([str(c.get_text()).split('_')[0] for c in a.get_xticklabels()])
        a.set_yticklabels([str(c.get_text()).split('_')[0] for c in a.get_yticklabels()])

    fig.savefig(argv.sample_corrs, dpi=300, bbox_inches='tight')


def plot_clone_corrs(cn, rt, clone_cn_cols, clone_rt_cols, argv):
    ''' Create figure that shows the correlation between the clone CN and RT profiles across all samples. '''
    fig, ax = plt.subplots(1, 2, figsize=(12,6), tight_layout=True)

    # pairwise distance between copy number profiles
    cn_dist = pd.DataFrame(
        squareform(pdist(cn.T.loc[clone_cn_cols])),
        columns = clone_cn_cols,
        index = clone_cn_cols
    )
    # subtract the max value from the distance matrix to make the heatmap more intuitive
    # this way the larger the value, the more similar the profiles
    cn_dist = cn_dist.max().max() - cn_dist
    # mask the upper triangle of the heatmap
    mask = np.zeros_like(cn_dist, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    # plot the heatmap
    sns.heatmap(cn_dist, square=False, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, ax=ax[0])
    ax[0].set_title('Similarity in clone CN')

    # correlation between RT profiles
    rt_corrs = rt[clone_rt_cols].corr()
    mask = np.zeros_like(rt_corrs, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(rt_corrs, square=False, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, ax=ax[1])
    ax[1].set_title('Correlation in clone RT')

    # only include the sample_id prefix and clone_id letter in the xticklabels and yticklabels
    ax[0].set_xticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2]) for c in ax[0].get_xticklabels()])
    ax[0].set_yticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2]) for c in ax[0].get_yticklabels()])

    # naming convention for the clone RT profiles is different than the clone CN profiles
    ax[1].set_xticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2].replace('clone', '')) for c in ax[1].get_xticklabels()])
    ax[1].set_yticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[2].replace('clone', '')) for c in ax[1].get_yticklabels()])
 
    fig.savefig(argv.clone_corrs, bbox_inches='tight', dpi=300)


def plot_rx_corrs(rt, clone_rt_cols, argv):
    ''' Plot the correlation between treated and untreated RT profiles. '''
    fig, ax = plt.subplots(1, 1, figsize=(4,4), tight_layout=True)

    # correlation between RT profiles
    rt_corrs = rt[clone_rt_cols].corr()
    mask = np.zeros_like(rt_corrs, dtype=bool)
    mask[np.triu_indices_from(mask)] = True
    sns.heatmap(rt_corrs, square=False, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, ax=ax)
    ax.set_title('Correlation in U vs T pseudobulk RT')

    # only include the sample ID prefix and the T or U suffix in the xticklabels and yticklabels
    ax.set_xticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[-1]) for c in ax.get_xticklabels()])
    ax.set_yticklabels(['{} {}'.format(str(c.get_text()).split('_')[0], str(c.get_text()).split('_')[-1]) for c in ax.get_yticklabels()])

    fig.savefig(argv.rx_corrs, bbox_inches='tight', dpi=300)



def main():
    argv = get_args()

    # read in the pseudobulk RT profiles for all samples
    rt = read_rt_data(argv)

    # read in the pseudobulk CN profiles for all samples
    cn = read_cn_data(argv)

    # load the number of cells belonging to each sample & clone
    counts = pd.read_csv(argv.counts, sep='\t')

    # create a list of the sample cn profiles
    sample_cn_cols = [c for c in cn.columns if c.endswith('sample_cn')]
    # create a list of the sample rt profiles
    sample_rt_cols = [c for c in rt.columns if c.endswith('model_rep_state') and 'clone' not in c]
    # create a list of the clone cn profiles
    temp_clone_cn_cols = [c for c in cn.columns if c.endswith('cn') and 'sample' not in c]
    # create a list of the clone rt profiles
    temp_clone_rt_cols = [c for c in rt.columns if c.endswith('model_rep_state') and 'clone' in c]
    # create a list of the treated vs untreated rt profiles
    rx_rt_cols = [c for c in rt.columns if c.endswith('model_rep_state_T') or c.endswith('model_rep_state_U')]

    clone_cn_cols = []
    clone_rt_cols = []

    # loop through all clones and add to the column list if the number of cells is greater than 10
    for col in temp_clone_cn_cols:
        # extract the dataset and clone_id from the column name
        c = col.split('_')[2]
        d = col.split('_')[0]
        # find the number of cells corresponding to this dataset and clone
        print('col', col, 'clone', c, 'dataset', d)
        # skip SA1035 clone F as there is not count data for this clone
        if c == 'F' and d == 'SA1035':
            continue
        n = int(counts.query("dataset=='{}'".format(d)).query("clone_id=='{}'".format(c))['num_cells_s'].values[0])
        if n > 10:
            clone_cn_cols.append(col)
    
    # repeat this for clone rt columns
    for col in temp_clone_rt_cols:
        c = col.split('_')[2].replace('clone', '')
        d = col.split('_')[0]
        print('col', col, 'clone', c, 'dataset', d)
        # skip SA1035 clone F as there is not count data for this clone
        if c == 'F' and d == 'SA1035':
            continue
        n = int(counts.query("dataset=='{}'".format(d)).query("clone_id=='{}'".format(c))['num_cells_s'].values[0])
        if n > 10:
            clone_rt_cols.append(col)

    # make sure the same number of clones are being plotted for CN and RT
    assert(len(clone_cn_cols) == len(clone_rt_cols))

    # plot the correlation between the sample CN and RT profiles
    plot_sample_corrs(cn, rt, sample_cn_cols, sample_rt_cols, argv)

    # plot the correlation between the clone CN and RT profiles
    plot_clone_corrs(cn, rt, clone_cn_cols, clone_rt_cols, argv)

    # plot the correlation between treated and untreated RT profiles
    plot_rx_corrs(rt, rx_rt_cols, argv)


if __name__ == '__main__':
    main()
