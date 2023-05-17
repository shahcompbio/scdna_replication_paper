import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_cna_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('rt', type=str, help='table of RT pseudobulk profiles')
    p.add_argument('cn', type=str, help='table of CN pseudobulk profiles')
    p.add_argument('counts', type=str, help='table of cell counts per clone')
    p.add_argument('rep_col', type=str, help='column for replication state (relevant for RT pseudobulk profile column names)')
    p.add_argument('dataset', type=str, help='name of this dataset')
    p.add_argument('out_table', type=str, help='Table of the number of cells per cell cycle phase and clone')
    p.add_argument('out_png', type=str, help='violin plot of subclonal RT differences split by CNA type')

    return p.parse_args()


def stratify_loci(df, clones, ploidies):
    """ Split all loci based on whether they contain no CNAs, clonal CNAs, or subclonal CNAs. """
    # mark all the indices in df that belong to each of the 3 categories
    no_event = []
    clonal_event = []
    subclonal_event = []
    for i, row in df.iterrows():
        # look at the relative difference between this site overall ploidy
        cn_diff = row[clones] - ploidies
        x = cn_diff.values.tolist()[0]
        # no event when all values are 0
        if np.all(cn_diff==0):
            no_event.append(i)
        # clonal event when all clones have the same nonzero
        # difference from the ploidy
        elif x.count(x[0])==len(x):
            clonal_event.append(i)
        # subclonal event when clones have different variations
        # from their respective ploidies
        else:
            subclonal_event.append(i)

    # split input df into 3 separate dfs based on CNA clonality of each locus
    df1 = df.iloc[no_event]
    df2 = df.iloc[clonal_event]
    df3 = df.iloc[subclonal_event]

    print(df1.shape, df2.shape, df3.shape)

    return df1, df2, df3


def subclonal_rt_diffs(df3, cols, clones, ploidies, argv):
    """ Find the RT difference between clones with a subclonal CNA and the reference (CN unaltered) clones at the same locus. """
    cna_rt_diffs = []
    for i, row in df3.iterrows():
        cn_diff = row[clones] - ploidies
        x = cn_diff.values.tolist()[0]
        unaltered_cols = []
        gain_cols = []
        loss_cols = []
        
        # get rt columns that correspond to clones with gain, loss, or no CNA at this locus
        for c, clone_id in enumerate(clones):
            if x[c]==0:
                unaltered_cols.append(cols[c])
            elif x[c] < 0:
                loss_cols.append(cols[c])
            else:
                gain_cols.append(cols[c])
       
        # define the reference RT for clones that don't have a CNA at this locus
        if len(unaltered_cols) > 0:
            avg_unaltered_rt = np.mean(row[unaltered_cols])
        else:
            # use the sample pseudbulk when there are subclonal CNAs in all clones
            avg_unaltered_rt = row['pseudobulk_{}'.format(argv.rep_col)]
        
        if len(gain_cols) > 0:
            avg_gain_rt = np.mean(row[gain_cols])
            gain_rt_diff = avg_gain_rt - avg_unaltered_rt
            temp_df = pd.DataFrame({
                'chr': [row['chr']], 'start': [row['start']], 'end': [row['end']],
                'reference_rt': [avg_unaltered_rt], 'clone_rt': [avg_gain_rt],
                'clone_rt_diff': [gain_rt_diff], 'clone_cna_type': ['gain'],
            })
            cna_rt_diffs.append(temp_df)

        if len(loss_cols) > 0:
            avg_loss_rt = np.mean(row[loss_cols])
            loss_rt_diff = avg_loss_rt - avg_unaltered_rt
            temp_df = pd.DataFrame({
                'chr': [row['chr']], 'start': [row['start']], 'end': [row['end']],
                'reference_rt': [avg_unaltered_rt], 'clone_rt': [avg_loss_rt],
                'clone_rt_diff': [loss_rt_diff], 'clone_cna_type': ['loss'],
            })
            cna_rt_diffs.append(temp_df)
        

    cna_rt_diffs = pd.concat(cna_rt_diffs, ignore_index=True)

    return cna_rt_diffs


def no_cna_rt_diffs(df1, cols, clones, ploidies, argv):
    """ Find the difference between a random reference clone and all other clones at a locus with no CNAs. This will serve as a background distribution. """
    unaltered_rt_diffs = []
    for i, row in df1.iterrows():
        # use the first clone as a reference
        reference_col = cols[0]
        reference_rt = row[reference_col]
        for c, clone_id in enumerate(clones):
            if c > 0:
                # compute the RT difference between this clone and the ref clone
                clone_rt = row[cols[c]]
                rt_diff = clone_rt - reference_rt
                
                temp_df = pd.DataFrame({
                    'chr': [row['chr']], 'start': [row['start']], 'end': [row['end']],
                    'reference_rt': [reference_rt], 'clone_rt': [clone_rt],
                    'clone_rt_diff': [rt_diff], 'clone_cna_type': ['unaltered'],
                })
                unaltered_rt_diffs.append(temp_df)

    unaltered_rt_diffs = pd.concat(unaltered_rt_diffs, ignore_index=True)

    return unaltered_rt_diffs


def violins_with_pvals(df, x, y, hue, ax, box_pairs, order=None, test='t-test_ind', text_format='star', loc='inside', verbose=0, palette=None):
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, order=order, palette=palette)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test, order=order,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_clone_rt_diff_vs_cna_types(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the distribution of clone RT differences against the CNA type of that particular locus. '''
    cna_cmap = get_cna_cmap()
    x = "clone_cna_type"
    y = "clone_rt_diff"
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


def main():
    argv = get_args()
    cn = pd.read_csv(argv.cn)
    counts = pd.read_csv(argv.counts)

    # drop the sample and dataset level cn pseudobulks
    good_cn_cols = [c for c in cn.columns if c.startswith('clone') or c in ['chr', 'start', 'end']]
    cn = cn[good_cn_cols]

    # set chr column to category
    cn.chr = cn.chr.astype(str)
    cn.chr = cn.chr.astype('category')

    rt = pd.read_csv(argv.rt)

    # set chr column to category
    rt.chr = rt.chr.astype(str)
    rt.chr = rt.chr.astype('category')

    df = pd.merge(cn, rt)
    df['end'] = df['start'] + 500000 - 1

    # find all the columns containing clone RT profiles
    temp_clone_rt_cols = [x for x in df.columns if x.startswith('pseudobulk_clone') and x.endswith(argv.rep_col)]

    # repeat this for clone rt columns
    cols = []
    d = argv.dataset
    for col in temp_clone_rt_cols:
        c = col.split('_')[1].replace('clone', '')
        print('col', col, 'clone', c, 'dataset', d)
        try:
            n = int(counts.query("dataset=='{}'".format(d)).query("clone_id=='{}'".format(c))['num_cells_s'].values[0])
        except:
            print('cannot find clone')
            n = 0
        if n > 10:
            cols.append(col)

    # find all the columns containing clone CN profiles
    clones = [x.split('_')[1].replace('clone', 'clone_') for x in cols]

    # compute the CN ploidy of each clone
    ploidies = df[clones].mode()
    
    # split loci based on CNA clonality
    df1, df2, df3 = stratify_loci(df, clones, ploidies)

    # compute clone RT differences in loci with subclonal CNAs
    cna_rt_diffs = subclonal_rt_diffs(df3, cols, clones, ploidies, argv)

    # find background distribition of clone RT differences in loci with no CNAs
    unaltered_rt_diffs = no_cna_rt_diffs(df1, cols, clones, ploidies, argv)

    # concatenate into one dataframe
    cna_rt_diffs = pd.concat([cna_rt_diffs, unaltered_rt_diffs], ignore_index=True)

    # create violinplot with t-tests for significance
    fig, ax = plt.subplots(1, 1, figsize=(6,4))
    ax = plot_clone_rt_diff_vs_cna_types(cna_rt_diffs, ax)
    ax.set_title('RT shifts at subclonal CNAs - {}'.format(argv.dataset))
    ax.set_xlabel('Subclonall CNA type')
    ax.set_ylabel('Clone RT relative to reference\n<--later | earlier-->')
    fig.savefig(argv.out_png, bbox_inches='tight', dpi=300)

    # add dataset name to table before saving
    cna_rt_diffs['dataset'] = argv.dataset

    # save a table of the RT diffs at each subclonal CNA site for later analysis
    cna_rt_diffs.to_csv(argv.out_table, index=False)



if __name__ == '__main__':
    main()
