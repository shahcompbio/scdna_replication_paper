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
    p.add_argument('-b0', '--beta0', type=float, nargs='+', help='beta0 for each dataset')
    p.add_argument('-b1', '--beta1', type=float, nargs='+', help='beta1 for each dataset')
    p.add_argument('-ns', '--num_s', type=int, nargs='+', help='number of S-phase cells in each dataset')
    p.add_argument('-krc', '--kronos_rep_col', type=str, help='column containing the kronos model replication states')
    p.add_argument('-prc', '--pert_rep_col', type=str, help='column containing the pert model replication states')
    p.add_argument('-pcn', '--pert_cn_col', type=str, help='column containing the pert model copy number states')
    p.add_argument('-trc', '--true_rep_col', type=str, help='column containing the true replication states')
    p.add_argument('-tcn', '--true_cn_col', type=str, help='column containing the true copy number states')
    p.add_argument('-t', '--table', help='table containing all the cn and rep accuracies for each simulated dataset and model')

    return p.parse_args()


def compute_accuracies(df, 
                       true_rep_col='true_rep',
                       model_rep_col='model_rep_state',
                       true_cn_col='true_G1_state',
                       model_cn_col='model_cn_state'):
    # compute the per-bin accuracy of true vs inferred copy number and replication states
    rep_acc  = 1.0 - (sum(abs(df[true_rep_col] - df[model_rep_col])) / df.shape[0])
    if model_cn_col is not None:
        cn_acc  = 1.0 - (sum(abs(df[true_cn_col] - df[model_cn_col])) / df.shape[0])
    else:
        cn_acc = None  # return None for models that don't infer CN
    return rep_acc, cn_acc


def load_data(chunk):
    # load the S-phase cells from each model version for a given dataset chunk
    df1 = pd.read_csv(chunk['kronos_path'].values[0], sep='\t')
    df2 = pd.read_csv(chunk['clone_path'].values[0], sep='\t')
    df3 = pd.read_csv(chunk['comp_path'].values[0], sep='\t')
    return df1, df2, df3


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
        'beta1': argv.beta1,
        'num_s': argv.num_s
    })

    # load table with paths to model results
    input_df = pd.read_csv(argv.input, sep='\t')

    # merge paths into legend_df
    legend_df = pd.merge(legend_df, input_df)

    df = []
    i = 0
    for dataset, chunk in legend_df.groupby('dataset'):
        df1, df2, df3 = load_data(chunk)

        # compute cn and rep accuracy for each method
        kronos_rep_acc, kronos_cn_acc = compute_accuracies(df1, model_rep_col=argv.kronos_rep_col, model_cn_col=None, true_cn_col=argv.true_cn_col, true_rep_col=argv.true_rep_col)
        pyro_clone_rep_acc, pyro_clone_cn_acc = compute_accuracies(df2, model_rep_col=argv.pert_rep_col, model_cn_col=argv.pert_cn_col, true_cn_col=argv.true_cn_col, true_rep_col=argv.true_rep_col)
        pyro_comp_rep_acc, pyro_comp_cn_acc = compute_accuracies(df3, model_rep_col=argv.pert_rep_col, model_cn_col=argv.pert_cn_col, true_cn_col=argv.true_cn_col, true_rep_col=argv.true_rep_col)

        datatag = dataset.split('.')[0]

        # create a dataframe with the accuracies for this simulated dataset
        methods = ['Kronos', 'PERT clone', 'PERT comp.']
        rep_accs = [kronos_rep_acc, pyro_clone_rep_acc, pyro_comp_rep_acc]
        cn_accs = [kronos_cn_acc, pyro_clone_cn_acc, pyro_comp_cn_acc]
        temp_df = pd.DataFrame({
            'dataset': [dataset]*3, 'datatag': [datatag]*3,
            'alpha': [chunk['alpha'].values[0]]*3, 'lambda': [chunk['lambda'].values[0]]*3,
            'beta0': [chunk['beta0'].values[0]]*3, 'beta1': [chunk['beta1'].values[0]]*3, 'num_s': [chunk['num_s'].values[0]]*3,
            'cell_cna_rate': [chunk.cell_cna_rate.values[0]]*3, 'num_clones': [chunk.num_clones.values[0]]*3,
            'method': methods, 'rep_accuracy': rep_accs, 'cn_accuracy': cn_accs
        })

        # store this dataframe of results and move onto the next simulated dataset
        df.append(temp_df)
        i += 1

    df = pd.concat(df, ignore_index=True)

    # save the raw table for posterity
    df.to_csv(argv.table, sep='\t', index=False)


if __name__ == '__main__':
    main()
