from argparse import ArgumentParser
import pandas as pd


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, help='table with path to model results for all datasets and inference methods')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-A', '--A', type=float, nargs='+', help='A params for each dataset')
    p.add_argument('-cna', '--cell_cna_rate', type=float, nargs='+', help='cell cna prob for each dataset')
    p.add_argument('-nc', '--num_clones', type=int, nargs='+', help='number of clones for each dataset')
    p.add_argument('-l', '--lamb', type=float, nargs='+', help='negative binomial event probs lambda for each dataset')
    p.add_argument('-b0', '--beta0', type=float, nargs='+', help='0th beta term for GC bias polynomial in each dataset')
    p.add_argument('-b1', '--beta1', type=float, nargs='+', help='1st beta term for GC bias polynomial in each dataset')
    p.add_argument('-tp', '--true_phase_col', type=str, help='column containing the true cell cycle phase')
    p.add_argument('-pp', '--pert_phase_col', type=str, help='column containing the pert predicted cell cycle phase')
    p.add_argument('-lp', '--laks_phase_col', type=str, help='column containing the laks predicted cell cycle phase')
    p.add_argument('-tf', '--true_frac_rep', help='true fraction of replicated bins per cell')
    p.add_argument('-pf', '--pert_frac_rep', help='PERT inferred fraction of replicated bins per cell')
    p.add_argument('-t', '--table', help='table containing all the cn and rep accuracies for each simulated dataset and model')

    return p.parse_args()


def load_data(chunk):
    # load the cn data for cells predicted to be in s-, g1-, and lq-phases
    cn_s = pd.read_csv(chunk['s_path'].values[0], sep='\t')
    cn_g = pd.read_csv(chunk['g_path'].values[0], sep='\t')
    cn_lq = pd.read_csv(chunk['lq_path'].values[0], sep='\t')
    return cn_s, cn_g, cn_lq


def compute_phase_accuracy(df, true_col, pred_col, true_frac_rep=None, pert_frac_rep=None):
    ''' Loop through every row in df, and count the number of times the true and predicted phases match. '''
    phase_acc = 0
    num_cells = df.shape[0]
    for i in range(num_cells):
        # count the phase as being correct if the true and predicted phases match
        if df.iloc[i][true_col] == df.iloc[i][pred_col]:
            phase_acc += 1
        # also count the phase as being correct if the cell is true S-phase
        # and both the true and predcted fraction of replicated bins are less than 0.05 or greater than 0.95
        elif df.iloc[i][true_col] == 'S' and true_frac_rep is not None and pert_frac_rep is not None:
            if df.iloc[i][true_frac_rep] < 0.05 and df.iloc[i][pert_frac_rep] < 0.05:
                phase_acc += 1
            elif df.iloc[i][true_frac_rep] > 0.95 and df.iloc[i][pert_frac_rep] > 0.95:
                phase_acc += 1
    phase_acc /= num_cells
    return phase_acc


def main():
    argv = get_args()

    # build table that matches config file for each permuted dataset
    legend_df = pd.DataFrame({
        'dataset': argv.datasets,
        'alpha': argv.A,
        'cell_cna_rate': argv.cell_cna_rate,
        'num_clones': argv.num_clones,
        'lambda': argv.lamb,
        'beta0': argv.beta0,
        'beta1': argv.beta1
    })

    # load table with paths to model results
    input_df = pd.read_csv(argv.input, sep='\t')

    # merge paths into legend_df
    legend_df = pd.merge(legend_df, input_df)

    df = []
    i = 0
    for dataset, chunk in legend_df.groupby('dataset'):
        cn_s, cn_g, cn_lq = load_data(chunk)

        # subset all three dataframes to just 'cell_id', 'true_phase', and 'PERT_phase
        cn_s = cn_s[['cell_id', argv.true_phase_col, argv.pert_phase_col, argv.true_frac_rep, argv.pert_frac_rep, argv.laks_phase_col]].drop_duplicates()
        cn_g = cn_g[['cell_id', argv.true_phase_col, argv.pert_phase_col, argv.true_frac_rep, argv.pert_frac_rep, argv.laks_phase_col]].drop_duplicates()
        cn_lq = cn_lq[['cell_id', argv.true_phase_col, argv.pert_phase_col, argv.true_frac_rep, argv.pert_frac_rep, argv.laks_phase_col]].drop_duplicates()

        # merge the three dataframes into one
        temp_df = pd.concat([cn_s, cn_g, cn_lq])

        # count the fraction of cells in which the true and predicted phases match
        pert_phase_acc = compute_phase_accuracy(temp_df, argv.true_phase_col, argv.pert_phase_col, argv.true_frac_rep, argv.pert_frac_rep)
        laks_phase_acc = compute_phase_accuracy(temp_df, argv.true_phase_col, argv.laks_phase_col)

        # add columns denoting dataset simulation parameters
        datatag = dataset.split('.')[0]
        temp_df['datatag'] = datatag
        temp_df['dataset'] = dataset
        temp_df['alpha'] = chunk['alpha'].values[0]
        temp_df['lambda'] = chunk['lambda'].values[0]
        temp_df['beta0'] = chunk['beta0'].values[0]
        temp_df['beta1'] = chunk['beta1'].values[0]
        temp_df['cell_cna_rate'] = chunk.cell_cna_rate.values[0]
        temp_df['num_clones'] = chunk.num_clones.values[0]
        temp_df['PERT_phase_acc'] = pert_phase_acc
        temp_df['laks_phase_acc'] = laks_phase_acc

        # store this dataframe of results and move onto the next simulated dataset
        df.append(temp_df)
        i += 1

    df = pd.concat(df, ignore_index=True)

    # save the raw table for posterity
    df.to_csv(argv.table, sep='\t', index=False)


if __name__ == '__main__':
    main()
