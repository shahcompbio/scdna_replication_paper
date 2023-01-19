from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, help='table with path to model results for all datasets and inference methods')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-A', '--A', type=float, nargs='+', help='A params for each dataset')
    p.add_argument('-cna', '--cell_cna_rate', type=float, nargs='+', help='cell cna prob for each dataset')
    p.add_argument('-nc', '--num_clones', type=int, nargs='+', help='number of clones for each dataset')
    p.add_argument('-l', '--lamb', type=float, nargs='+', help='negative binomial event probs lambda for each dataset')
    p.add_argument('-brc', '--bulk_rep_col', type=str, help='column containing the bulk model replication states')
    p.add_argument('-prc', '--pyro_rep_col', type=str, help='column containing the pyro model replication states')
    p.add_argument('-pcn', '--pyro_cn_col', type=str, help='column containing the pyro model copy number states')
    p.add_argument('-trc', '--true_rep_col', type=str, help='column containing the true replication states')
    p.add_argument('-tcn', '--true_cn_col', type=str, help='column containing the true copy number states')
    p.add_argument('-t', '--table', help='table containing all the cn and rep accuracies for each simulated dataset and model')
    p.add_argument('-p1', '--plot1', help='plot showing all the cn and rep accuracies for each simulated dataset and model')
    p.add_argument('-p2', '--plot2', help='plot of parameter sweep across num_clones and cell_cna_rate')

    return p.parse_args()


def make_figure1(df, argv):
    """ First figure is a mixture of barplots, scatterplots, and violinplots showing the accuracies of the different methods. """
    fig, ax = plt.subplots(2, 4, figsize=(16, 8), tight_layout=True)

    # merge together the two supblots in the top left corner
    gs = ax[0, 0].get_gridspec()
    for a in ax[0, :2]:
        a.remove()
    axbig_top_row = fig.add_subplot(gs[0, 0:2])

    # merge together the two supblots in the bottom left corner
    gs = ax[1, 0].get_gridspec()
    for a in ax[1, :2]:
        a.remove()
    axbig_bottom_row = fig.add_subplot(gs[1, 0:2])

    # showing rep accuracies on the top row
    # barplots and scatterplots of cn and rep accuracies for each method, across different simulation params
    sns.barplot(data=df, x='datatag', y='rep_accuracy', hue='method', ax=axbig_top_row)

    # scatterplots which use A and lambda as the size params
    sns.scatterplot(data=df, x='cell_cna_rate', y='rep_accuracy', hue='method', size='alpha', style='num_clones', ax=ax[0, 2])
    sns.scatterplot(data=df, x='cell_cna_rate', y='rep_accuracy', hue='method', size='lambda', style='num_clones', ax=ax[0, 3])

    # showing cn accuracies on the bottom row
    sns.barplot(data=df, x='datatag', y='cn_accuracy', hue='method', ax=axbig_bottom_row)
    
    sns.scatterplot(data=df, x='cell_cna_rate', y='cn_accuracy', hue='method', size='alpha', style='num_clones', ax=ax[1, 2])
    sns.scatterplot(data=df, x='cell_cna_rate', y='cn_accuracy', hue='method', size='lambda', style='num_clones', ax=ax[1, 3])

    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def violins_with_pvals(df, x, y, hue, ax, box_pairs, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    """ Create a violinplot with p-values annotated. """
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test,
                        text_format=text_format, loc=loc, verbose=verbose)


def plot_cna_rate_rep_acc(df, ax, n=1, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the replications state accuracy vs cna rate at a fixed number of clones (n). '''
    x = "cell_cna_rate"
    y = "rep_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('num_clones=={}'.format(n))
    box_pairs = [
        ((0.05, "PERT clone"), (0.05, "PERT comp.")),
        ((0.02, "PERT clone"), (0.02, "PERT comp.")),
        ((0.02, "PERT comp."), (0.00, "PERT comp.")),
        ((0.02, "PERT comp."), (0.05, "PERT comp.")),
        ((0.00, "PERT comp."), (0.05, "PERT comp.")),
        ((0.05, "Dileep"), (0.05, "PERT comp.")),
        ((0.05, "PERT clone"), (0.05, "Dileep")),
        ((0.02, "Dileep"), (0.02, "PERT comp.")),
        ((0.02, "PERT clone"), (0.02, "Dileep")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test, text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('{} clone(s)'.format(n))
    return ax


def plot_cna_rate_cn_acc(df, ax, n=1, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the cn accuracy vs cna rate at a fixed number of clones (n). '''
    x = "cell_cna_rate"
    y = "cn_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('num_clones=={}'.format(n))
    box_pairs = [
        ((0.05, "PERT clone"), (0.05, "PERT comp.")),
        ((0.02, "PERT clone"), (0.02, "PERT comp.")),
        ((0.02, "PERT comp."), (0.00, "PERT comp.")),
        ((0.02, "PERT comp."), (0.05, "PERT comp.")),
        ((0.00, "PERT comp."), (0.05, "PERT comp.")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('{} clone(s)'.format(n))
    return ax


def plot_clone_effect_rep_acc(df, ax, rate=0.02, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the replication accuracy vs number of clones at a fixed cell cna rate. '''
    x = "num_clones"
    y = "rep_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('cell_cna_rate=={}'.format(rate)).query('num_clones<4')
    box_pairs = [
        ((1, "Dileep"), (3, "Dileep")),
        ((1, "PERT comp."), (3, "PERT comp.")),
        ((3, "Dileep"), (3, "PERT comp.")),
        ((3, "PERT clone"), (3, "Dileep")),
        ((3, "PERT clone"), (3, "PERT comp.")),
        ((1, "Dileep"), (1, "PERT comp.")),
        ((1, "PERT clone"), (1, "Dileep")),
        ((1, "PERT clone"), (1, "PERT comp.")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('Cell CNA rate {}'.format(rate))
    return ax


def plot_clone_effect_cn_acc(df, ax, rate=0.02, test='t-test_ind', text_format='star', loc='inside', verbose=0):
    ''' Plot the CN accuracy vs number of clones at a fixed cell cna rate. '''
    x = "num_clones"
    y = "cn_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.7').query('cell_cna_rate=={}'.format(rate)).query('num_clones<4')
    box_pairs = [
        ((1, "PERT comp."), (3, "PERT comp.")),
        ((3, "PERT clone"), (3, "PERT comp.")),
        ((1, "PERT clone"), (1, "PERT comp.")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose)
    ax.set_title('Cell CNA rate {}'.format(rate))
    return ax


def make_figure2(df, argv):
    """ Figure showing parameter sweep across num_clones and cell_cna_rate """
    fig, ax = plt.subplots(2, 5, figsize=(20, 8), tight_layout=True)
    ax = ax.flatten()

    # top row shows accuracy at predicting replication states
    # show effect of varying cna rate at fixed number of clones
    plot_cna_rate_rep_acc(df, ax[0], n=1)
    plot_cna_rate_rep_acc(df, ax[1], n=3)
    # show effect of varying number of clones at fixed cna rates
    plot_clone_effect_rep_acc(df, ax[2], rate=0.0)
    plot_clone_effect_rep_acc(df, ax[3], rate=0.02)
    plot_clone_effect_rep_acc(df, ax[4], rate=0.05)

    # bottom row shows accuracy at predicting CN states
    # show effect of varying cna rate at fixed number of clones
    plot_cna_rate_cn_acc(df, ax[5], n=1)
    plot_cna_rate_cn_acc(df, ax[6], n=3)
    # show effect of varying number of clones at fixed cna rates
    plot_clone_effect_cn_acc(df, ax[7], rate=0.0)
    plot_clone_effect_cn_acc(df, ax[8], rate=0.02)
    plot_clone_effect_cn_acc(df, ax[9], rate=0.05)

    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)



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
        'alpha': argv.A,
        'cell_cna_rate': argv.cell_cna_rate,
        'num_clones': argv.num_clones,
        'lambda': argv.lamb
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
        methods = ['Dileep', 'PERT clone', 'PERT comp.']
        rep_accs = [bulk_rep_acc, pyro_clone_rep_acc, pyro_comp_rep_acc]
        cn_accs = [bulk_cn_acc, pyro_clone_cn_acc, pyro_comp_cn_acc]
        temp_df = pd.DataFrame({
            'dataset': [dataset]*3, 'datatag': [datatag]*3,
            'alpha': [chunk['alpha'].values[0]]*3, 'lambda': [chunk['lambda'].values[0]]*3,
            'cell_cna_rate': [chunk.cell_cna_rate.values[0]]*3, 'num_clones': [chunk.num_clones.values[0]]*3,
            'method': methods, 'rep_accuracy': rep_accs, 'cn_accuracy': cn_accs
        })

        # store this dataframe of results and move onto the next simulated dataset
        df.append(temp_df)
        i += 1

    df = pd.concat(df, ignore_index=True)

    # make the first figure
    make_figure1(df, argv)

    # rename lambda column to avoid conflict with python keyword
    df.rename(columns={'lambda': 'lamb'}, inplace=True)

    # make the second figure
    make_figure2(df, argv)

    # save the raw table for posterity
    df.to_csv(argv.table, sep='\t', index=False)


if __name__ == '__main__':
    main()
