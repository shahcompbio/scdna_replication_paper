from argparse import ArgumentParser
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, help='table with path to model results for all datasets and inference methods')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-A', '--A', type=float, nargs='+', help='A params for each dataset')
    p.add_argument('-cna', '--cell_cna_prob', type=float, nargs='+', help='cell cna prob for each dataset')
    p.add_argument('-nc', '--num_clones', type=int, nargs='+', help='number of clones for each dataset')
    p.add_argument('-nbr', '--nb_r', type=float, nargs='+', help='negative binomial rate for each dataset')
    p.add_argument('-brc', '--bulk_rep_col', type=str, help='column containing the bulk model replication states')
    p.add_argument('-prc', '--pyro_rep_col', type=str, help='column containing the pyro model replication states')
    p.add_argument('-pcn', '--pyro_cn_col', type=str, help='column containing the pyro model copy number states')
    p.add_argument('-trc', '--true_rep_col', type=str, help='column containing the true replication states')
    p.add_argument('-tcn', '--true_cn_col', type=str, help='column containing the true copy number states')
    p.add_argument('-t', '--table', help='table containing all the cn and rep accuracies for each simulated dataset and model')
    p.add_argument('-p', '--plot', help='plot showing all the cn and rep accuracies for each simulated dataset and model')

    return p.parse_args()


def make_figure(df, argv):
    fig, ax = plt.subplots(2, 5, figsize=(25, 8), tight_layout=True)
    ax = ax.flatten()

    # barplots and scatterplots of cn and rep accuracies for each model, across different simulation params
    # showing rep accuracies on the top row
    sns.barplot(data=df, x='datatag', y='rep_accuracy', hue='model', ax=ax[0])
    # ax[0].set_ylim(0.7, 1.02)

    sns.violinplot(data=df.query("num_clones==1"), x='cell_cna_prob', y='rep_accuracy', hue='model', ax=ax[1])
    ax[1].set_title('Diploid')
    # ax[1].set_ylim(0.7, 1.02)
    
    sns.violinplot(data=df.query("cell_cna_prob==0.02"), x='num_clones', y='rep_accuracy', hue='model', ax=ax[2])
    ax[2].set_title('Cell CNA rate 2%')
    # ax[2].set_ylim(0.7, 1.02)

    sns.scatterplot(data=df, x='cell_cna_prob', y='rep_accuracy', hue='model', size='A', style='num_clones', ax=ax[3])

    sns.scatterplot(data=df, x='cell_cna_prob', y='rep_accuracy', hue='model', size='nb_r', style='num_clones', ax=ax[4])

    # showing cn accuracies on the bottom row
    sns.barplot(data=df, x='datatag', y='cn_accuracy', hue='model', ax=ax[5])
    # ax[5].set_ylim(0.7, 1.02)

    sns.violinplot(data=df.query("num_clones==1"), x='cell_cna_prob', y='cn_accuracy', hue='model', ax=ax[6])
    ax[6].set_title('Diploid')
    # ax[6].set_ylim(0.7, 1.02)

    sns.violinplot(data=df.query("cell_cna_prob==0.02"), x='num_clones', y='cn_accuracy', hue='model', ax=ax[7])
    ax[7].set_title('Cell CNA rate 2%')
    # ax[7].set_ylim(0.7, 1.02)
    
    sns.scatterplot(data=df, x='cell_cna_prob', y='cn_accuracy', hue='model', size='A', style='num_clones', ax=ax[8])
    
    sns.scatterplot(data=df, x='cell_cna_prob', y='cn_accuracy', hue='model', size='nb_r', style='num_clones', ax=ax[9])

    fig.savefig(argv.plot, bbox_inches='tight')



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
    df1 = pd.read_csv(chunk['bulk_path'].values[0], sep='\t')
    df2 = pd.read_csv(chunk['clone_path'].values[0], sep='\t')
    df3 = pd.read_csv(chunk['comp_path'].values[0], sep='\t')
    return df1, df2, df3


def main():
    argv = get_args()

    # build table that matches config file for each permuted dataset
    legend_df = pd.DataFrame({
        'dataset': argv.datasets,
        'A': argv.A,
        'cell_cna_prob': argv.cell_cna_prob,
        'num_clones': argv.num_clones,
        'nb_r': argv.nb_r
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
        bulk_rep_acc, bulk_cn_acc = compute_accuracies(df1, model_rep_col=argv.bulk_rep_col, model_cn_col=None, true_cn_col=argv.true_cn_col, true_rep_col=argv.true_rep_col)
        pyro_clone_rep_acc, pyro_clone_cn_acc = compute_accuracies(df2, model_rep_col=argv.pyro_rep_col, model_cn_col=argv.pyro_cn_col, true_cn_col=argv.true_cn_col, true_rep_col=argv.true_rep_col)
        pyro_comp_rep_acc, pyro_comp_cn_acc = compute_accuracies(df3, model_rep_col=argv.pyro_rep_col, model_cn_col=argv.pyro_cn_col, true_cn_col=argv.true_cn_col, true_rep_col=argv.true_rep_col)

        datatag = dataset.split('.')[0]

        # create a dataframe with the accuracies for this simulated dataset
        models = ['bulk', 'pyro_clone', 'pyro_comp']
        rep_accs = [bulk_rep_acc, pyro_clone_rep_acc, pyro_comp_rep_acc]
        cn_accs = [bulk_cn_acc, pyro_clone_cn_acc, pyro_comp_cn_acc]
        temp_df = pd.DataFrame({
            'dataset': [dataset]*3, 'datatag': [datatag]*3,
            'A': [chunk.A.values[0]]*3, 'nb_r': [chunk.nb_r.values[0]]*3,
            'cell_cna_prob': [chunk.cell_cna_prob.values[0]]*3, 'num_clones': [chunk.num_clones.values[0]]*3,
            'model': models, 'rep_accuracy': rep_accs, 'cn_accuracy': cn_accs
        })

        # store this dataframe of results and move onto the next simulated dataset
        df.append(temp_df)
        i += 1

    df = pd.concat(df, ignore_index=True)

    # make the desired figure
    make_figure(df, argv)

    # save the raw table for posterity
    df.to_csv(argv.table, sep='\t', index=False)


if __name__ == '__main__':
    main()
