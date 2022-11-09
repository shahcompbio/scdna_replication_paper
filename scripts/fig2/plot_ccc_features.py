import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('pyro_input', type=str, help='pyro model output for s-phase cells')
    p.add_argument('bulk_input', type=str, help='bulk model output for s-phase cells')
    p.add_argument('g1_input', type=str, help='g1-phase cells with ccc features')
    p.add_argument('frac_rt_col', type=str, help='column name for the fraction of replicated bins in a cell')
    p.add_argument('pyro_rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('bulk_rep_col', type=str, help='column name for replicated status of each bin in bulk model')
    p.add_argument('plot1', type=str, help='figure showing histogram of ccc features based on cell cycle phase and extreme_cell_frac status')
    p.add_argument('plot2', type=str, help='figure showing scatterplot of ccc features based on cell cycle phase and extreme_cell_frac status')

    return p.parse_args()


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state'):
    ''' Compute the fraction of replicated bins for all cells in `cn`, noting which cells have extreme values. '''
    cn['extreme_cell_frac'] = False
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn[rep_state_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac

        if temp_frac > 0.95 or temp_frac < 0.05:
            cn.loc[cell_cn.index, 'extreme_cell_frac'] = True
    return cn


def compute_ecf(cn_pyro, cn_bulk, frac_rt_col='cell_frac_rep'):
    ''' Find which cells have extreme cell fraction of replicated bins for both the pyro and bulk models. '''
    cell_metric_cols = [
        'cell_id', frac_rt_col, 
        'extreme_cell_frac', 'breakpoints', 
        'num_reads', 'madn', 'lrs', 
        'corrected_madn', 'corrected_breakpoints'
    ]
    cell_metrics_pyro = cn_pyro[cell_metric_cols].drop_duplicates()
    cell_metrics_bulk = cn_bulk[cell_metric_cols].drop_duplicates()

    ecf_pyro = cell_metrics_pyro[['cell_id', 'extreme_cell_frac']].drop_duplicates().reset_index(drop=True).rename(columns={'extreme_cell_frac': 'ecf_pyro'})
    ecf_bulk = cell_metrics_bulk[['cell_id', 'extreme_cell_frac']].drop_duplicates().reset_index(drop=True).rename(columns={'extreme_cell_frac': 'ecf_bulk'})
    ecf = pd.merge(ecf_pyro, ecf_bulk)

    return ecf


def get_cell_metrics(cn_pyro, cn_g1, ecf, frac_rt_col='cell_frac_rep'):
    # merge information about where ecf cells are located
    cn_pyro = pd.merge(cn_pyro, ecf)

    # use ecf status to nominate "bad" S-phase cells that are likely nonreplicating
    cn_pyro['cell_cycle_state'] = 'S'
    cn_g1['cell_cycle_state'] = 'G1/2'
    print(ecf['ecf_pyro'].value_counts())
    print(cn_pyro['ecf_pyro'].value_counts())
    print(cn_pyro['ecf_bulk'].value_counts())
    extreme_cn = cn_pyro.loc[(cn_pyro['ecf_pyro']==True) | (cn_pyro['ecf_bulk']==True)]
    print(extreme_cn.shape)
    print(extreme_cn.head())
    cn_pyro.loc[extreme_cn.index, 'cell_cycle_state'] = 'bad-S'
    print(cn_pyro.cell_cycle_state.value_counts())

    # get cell-level metrics with the updated cell cycle labels
    cell_metric_cols = [
        'cell_id', 'cell_cycle_state', 'breakpoints', 
        'num_reads', 'madn', 'lrs', frac_rt_col, 
        'corrected_madn', 'corrected_breakpoints'
    ]
    cell_metrics_pyro = cn_pyro[cell_metric_cols].drop_duplicates()
    cell_metrics_g1 = cn_g1[cell_metric_cols].drop_duplicates()

    cell_metrics = pd.concat([cell_metrics_pyro, cell_metrics_g1], ignore_index=True)

    return cell_metrics


def plot_ccc_distributions(cell_metrics, argv):
    # plot histogram of each feature split by cell cycle state
    fig, ax = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    sns.histplot(data=cell_metrics, x='corrected_madn', hue='cell_cycle_state', multiple='stack', ax=ax[0])
    sns.histplot(data=cell_metrics, x='lrs', hue='cell_cycle_state', multiple='stack', ax=ax[1])
    sns.histplot(data=cell_metrics, x='corrected_breakpoints', hue='cell_cycle_state', multiple='stack', ax=ax[2])

    fig.savefig(argv.plot1, bbox_inches='tight')

    # plot pairwise scatterplots of all three features colored by cell cycle cstate
    fig, ax = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(data=cell_metrics, x='corrected_madn', y='lrs', hue='cell_cycle_state', alpha=0.3, ax=ax[0])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='lrs', hue='cell_cycle_state', alpha=0.3, ax=ax[1])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='corrected_madn', hue='cell_cycle_state', alpha=0.3, ax=ax[2])
    sns.scatterplot(data=cell_metrics, x='corrected_madn', y='lrs', hue=argv.frac_rt_col, alpha=0.3, ax=ax[3])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='lrs', hue=argv.frac_rt_col, alpha=0.3, ax=ax[4])
    sns.scatterplot(data=cell_metrics, x='corrected_breakpoints', y='corrected_madn', hue=argv.frac_rt_col, alpha=0.3, ax=ax[5])

    fig.savefig(argv.plot2, bbox_inches='tight')


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn_pyro = pd.read_csv(argv.pyro_input, sep='\t')
    cn_bulk = pd.read_csv(argv.bulk_input, sep='\t')
    cn_g1 = pd.read_csv(argv.g1_input, sep='\t')

    # compute the fraction of replicated bins within each cell
    cn_pyro = compute_cell_frac(cn_pyro, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.pyro_rep_col)
    cn_bulk = compute_cell_frac(cn_bulk, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.bulk_rep_col)

    # dummy column for G1-phase cells to have frac_rep of -1
    cn_g1[argv.frac_rt_col] = -1

    # label which cells have extreme_cell_frac for both the pyro and bulk
    # versions of inferring scRT states
    ecf = compute_ecf(cn_pyro, cn_bulk, frac_rt_col=argv.frac_rt_col)

    # use ecf status to come up with good vs bad S-phase labels in cell_cycle_state column
    cell_metrics = get_cell_metrics(cn_pyro, cn_g1, ecf, frac_rt_col=argv.frac_rt_col)

    # generate and save plots
    plot_ccc_distributions(cell_metrics, argv)
