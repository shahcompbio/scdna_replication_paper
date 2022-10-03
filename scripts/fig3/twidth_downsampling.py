from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scdna_replication_tools.calculate_twidth import calculate_twidth
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('-cn', '--cn_s', type=str, nargs='+', help='model results from pyro with composite cn prior')
    p.add_argument('-d', '--dataset', type=str, nargs='+')
    p.add_argument('-l', '--cell_type_labels', type=str, nargs='+')
    p.add_argument('-frc', '--frac_rt_col', help='column denoting the fraction of replicated loci per cell (its time in S-phase)')
    p.add_argument('-rc', '--rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('--table', help='table of all the computed t_width values')
    p.add_argument('--plot', help='T-width curve using inferred scRT')

    return p.parse_args()


def load_data(argv):
    # load long-form cn data for each dataset passed into the model
    i = 0
    df = []
    for path in argv.cn_s:
        temp_df = pd.read_csv(path, sep='\t')
        temp_df = temp_df[['chr', 'start', 'cell_id', 'clone_id', argv.frac_rt_col, argv.rep_col]]
        temp_df['dataset'] = argv.dataset[i]
        temp_df['cell_type'] = argv.cell_type_labels[i]
        df.append(temp_df)
        i += 1

    # concatenate into one df
    df = pd.concat(df, ignore_index=True)
    return df


def compute_twidth_de_novo(df, rep_col='model_rep_state', frac_rep_col='cell_frac_rep', dataset=None):
    # note the number of cells in this df
    num_cells = len(df.cell_id.unique())
    
    # compute the pseudobulk rt profile for all the cells in the df
    bulk_rt = compute_pseudobulk_rt_profiles(
        df, rep_col, output_col='pseduobulk', time_col='hours', clone_col='clone_id'
    )
    df = pd.merge(df, bulk_rt)
    
    # compute time from scheduled replication column
    df['time_from_scheduled_rt'] = df['pseduobulk_hours'] - (df[frac_rep_col] * 10.0)
    
    # dataframe storing computed T-width values
    t_width_df = []

    # compute and plot twidth for per_cell==False
    Tw, _, _, _, _, _ = calculate_twidth(
        df, tfs_col='time_from_scheduled_rt', rs_col=rep_col, per_cell=False, curve='sigmoid'
    )
    t_width_df.append(pd.DataFrame({'dataset': [dataset], 'num_cells': [num_cells], 'per_cell': [False], 'T-width': [Tw]}))

    # compute and plot twidth for per_cell==True
    Tw, _, _, _, _, _ = calculate_twidth(
        df, tfs_col='time_from_scheduled_rt', rs_col=rep_col, per_cell=True, curve='sigmoid'
    )
    t_width_df.append(pd.DataFrame({'dataset': [dataset], 'num_cells': [num_cells], 'per_cell': [True], 'T-width': [Tw]}))

    t_width_df = pd.concat(t_width_df, ignore_index=True)
    return t_width_df


def subset_cells(df, num_cells=None, frac_cells=None):
    # subset df to a given number (or fraction) of cells
    assert (num_cells is not None) or (frac_cells is not None)
    if num_cells is None:
        num_cells = frac_cells * len(df.cell_id.unique())
    
    good_cells = np.random.choice(df.cell_id.unique(), num_cells, replace=False)
    out_df = df.loc[df['cell_id'].isin(good_cells)]
    
    return out_df



def main():
    argv = get_args()
    
    # load all the s-phase cells from all datasets into one dataframe
    df = load_data(argv)

    # test population sizes of 25, 50, 75, ..., 375 cells
    num_cells = np.arange(25, 400, 25)
    t_width_df = []
    for dataset, chunk_df in df.groupby('dataset'):
        # only subsample the datasets with >375 S-phase cells
        if len(chunk_df.cell_id.unique()) > max(num_cells):
            # compute Twidth at smaller subsets
            for n in num_cells:
                for i in range(5):
                    print('computing twidth for {} at size {} cells and sample {}'.format(dataset, n, i))
                    temp_df = subset_cells(chunk_df, num_cells=n)  # subset this dataset to the desired number of cells
                    temp_tw_df = compute_twidth_de_novo(temp_df, rep_col=argv.rep_col, frac_rep_col=argv.frac_rt_col, dataset=dataset)
                    t_width_df.append(temp_tw_df)
        # compute T-width with no subsampling
        temp_tw_df = compute_twidth_de_novo(chunk_df, rep_col=argv.rep_col, frac_rep_col=argv.frac_rt_col, dataset=dataset)
        t_width_df.append(temp_tw_df)
    t_width_df = pd.concat(t_width_df, ignore_index=True)

    # save figure of twidth curves
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)
    # cell_type_order = ['WT', 'TP53-/-', 'TP53-/-, BRCA1+/-', 'TP53-/-, BRCA1-/-', 'TP53-/-, BRCA2+/-', 'TP53-/-, BRCA2-/-']
    
    sns.scatterplot(data=df.query('per_cell==False'), y='num_cells', x='T-width', hue='cell_type', ax=ax[0])
    ax[0].set_title('Downsampling T-width')

    sns.scatterplot(data=df.query('per_cell==False'), y='num_cells', x='T-width', hue='cell_type', ax=ax[1])
    ax[1].set_title('Downsampling T-width (per-cell)')

    fig.savefig(argv.plot, bbox_inches='tight')

    # save a table of all the computed T-width values
    t_width_df.to_csv(argv.table, sep='\t', index=False)


if __name__ == '__main__':
    main()
