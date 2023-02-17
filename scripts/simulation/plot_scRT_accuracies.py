from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='table containing all the cn and rep accuracies for each simulated dataset and model')
    p.add_argument('plot1', help='plot showing all the cn and rep accuracies for each simulated dataset and model')
    p.add_argument('plot2', help='plot of parameter sweep across num_clones and cell_cna_rate')

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
        ((0.05, "Kronos"), (0.05, "PERT comp.")),
        ((0.05, "PERT clone"), (0.05, "Kronos")),
        ((0.02, "Kronos"), (0.02, "PERT comp.")),
        ((0.02, "PERT clone"), (0.02, "Kronos")),
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
        ((1, "Kronos"), (3, "Kronos")),
        ((1, "PERT comp."), (3, "PERT comp.")),
        ((3, "Kronos"), (3, "PERT comp.")),
        ((3, "PERT clone"), (3, "Kronos")),
        ((3, "PERT clone"), (3, "PERT comp.")),
        ((1, "Kronos"), (1, "PERT comp.")),
        ((1, "PERT clone"), (1, "Kronos")),
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



def main():
    argv = get_args()

    # load the data
    df = pd.read_csv(argv.input, sep='\t')

    # make the first figure
    make_figure1(df, argv)

    # rename lambda column to avoid conflict with python keyword
    df.rename(columns={'lambda': 'lamb'}, inplace=True)

    # make the second figure
    make_figure2(df, argv)



if __name__ == '__main__':
    main()
