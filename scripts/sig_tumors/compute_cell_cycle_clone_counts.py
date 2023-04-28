import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', type=str, help='full df for S-phase cells')
    p.add_argument('cn_g', type=str, help='full df for G1/2-phase cells')
    p.add_argument('out_csv', type=str, help='Table of the number of cells per cell cycle phase and clone')

    return p.parse_args()


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
    

def library_wrapper_fracs_and_pvals(df):
    ''' Compute clone cell cycle fractions within each library '''
    df_out = []
    for library_id, chunk in df.groupby('library_id'):
        temp_out = compute_fracs_and_pvals(chunk)
        temp_out['library_id'] = library_id
        df_out.append(temp_out)
    df_out = pd.concat(df_out, ignore_index=True)
    return df_out


def clone_spf_analysis(cn_s, cn_g, argv):
    # add column to denote phase of each df
    cn_s['cell_cycle'] = 'S'
    cn_g['cell_cycle'] = 'G1/2'

    # # rename the clone_id column to match the cells in the tree
    # cn_s['clone_id'] = cn_s['assigned_clone_id']
    # cn_g['clone_id'] = cn_g['assigned_clone_id']

    # concatenate all cells into one df but only store relevant columns
    coi = ['cell_id', 'cell_cycle', 'clone_id', 'library_id']
    df = pd.concat([cn_s[coi], cn_g[coi]], ignore_index=True)
    df = df.sort_values(['clone_id', 'cell_cycle', 'library_id']).drop_duplicates(ignore_index=True)

    # count the number of cells belonging to each clone and cell cycle phase
    df2 = df[['cell_cycle', 'clone_id', 'library_id']].value_counts().to_frame(name='num_cells').reset_index().sort_values(['clone_id', 'cell_cycle']).reset_index(drop=True)
    
    # remove clones with "None" ID
    df2 = df2.loc[df2['clone_id']!='None']

    clones = df2['clone_id'].unique()
    libraries = df2['library_id'].unique()

    # add all the counts for clones+phase combos with no cells
    absent_df = []
    for library_id in libraries:
        for clone_id in clones:
            if clone_id not in df2.query("cell_cycle=='S'").query("library_id=='{}'".format(library_id))['clone_id'].values:
                absent_df.append(pd.DataFrame({
                    'cell_cycle': ['S'], 'clone_id': [clone_id], 'library_id': [library_id], 'num_cells': [0]
                }))
            if clone_id not in df2.query("cell_cycle=='G1/2'").query("library_id=='{}'".format(library_id))['clone_id'].values:
                absent_df.append(pd.DataFrame({
                    'cell_cycle': ['G1/2'], 'clone_id': [clone_id], 'library_id': [library_id], 'num_cells': [0]
                }))
    
    # concatenate into one dataframe
    if len(absent_df) > 0:
        absent_df = pd.concat(absent_df, ignore_index=True)
        df2 = pd.concat([df2, absent_df], ignore_index=True)

    # compute cell cycle fractions per clone & library along with p-values
    df2 = library_wrapper_fracs_and_pvals(df2)

    # Bonferroni correction of p-values by the total number of hypotheses tested
    # x2 because we're doing a two-sided t-test
    df2['positive_p_adj'] = df2['positive_p'] * 2 * df2.shape[0]
    df2['negative_p_adj'] = df2['negative_p'] * 2 * df2.shape[0]

    # save table used to generate spf figure
    df2.to_csv(argv.out_csv, index=False)


def main():
    argv = get_args()
    # load long-form dataframes from different cell cycle phases
    cn_s = pd.read_csv(argv.cn_s)
    cn_g = pd.read_csv(argv.cn_g)

    cn_s.chr = cn_s.chr.astype('str')
    cn_s.chr = cn_s.chr.astype('category')
    cn_g.chr = cn_g.chr.astype('str')
    cn_g.chr = cn_g.chr.astype('category')

    # plot S-phase fractions at the clone and sample levels
    # save both the plot and table used to make the plot
    clone_spf_analysis(cn_s, cn_g, argv)


if __name__ == '__main__':
    main()

