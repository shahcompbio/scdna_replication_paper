from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_methods_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='table containing all the cn and rep accuracies for each simulated dataset and model')
    p.add_argument('plot1', help='plot showing all the cn and rep accuracies for each simulated dataset and model')
    p.add_argument('plot2', help='plot of parameter sweep across num_clones and cell_cna_rate')

    return p.parse_args()


def violins_with_pvals(df, x, y, hue, ax, box_pairs, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    """ Create a violinplot with p-values annotated. """
    palette = get_methods_cmap()
    sns.violinplot(data=df, x=x, y=y, hue=hue, ax=ax, palette=palette, saturation=1, linewidth=1, hue_order=hue_order)
    add_stat_annotation(ax, data=df, x=x, y=y, hue=hue,
                        box_pairs=box_pairs, test=test,
                        text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)


def plot_cna_rate_rep_acc(df, ax, n=1, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the replications state accuracy vs cna rate at a fixed number of clones (n). '''
    x = "cell_cna_rate"
    y = "rep_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.75').query('num_clones=={}'.format(n)).query('beta0==1.2')
    box_pairs = [
        ((0.05, "PERT clone"), (0.05, "PERT comp.")),
        ((0.02, "PERT clone"), (0.02, "PERT comp.")),
        ((0.05, "Kronos"), (0.05, "PERT comp.")),
        ((0.05, "PERT clone"), (0.05, "Kronos")),
        ((0.02, "Kronos"), (0.02, "PERT comp.")),
        ((0.02, "PERT clone"), (0.02, "Kronos")),
        ((0.00, "PERT comp."), (0.00, "PERT clone")),
        ((0.00, "PERT comp."), (0.00, "Kronos")),
        ((0.00, "Kronos"), (0.00, "PERT clone")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test, text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across cell CNA rate\n# of clones={}'.format(n))
    ax.set_xlabel('Cell CNA rate')
    ax.set_ylabel('Replication state accuracy')
    return ax


def plot_cna_rate_cn_acc(df, ax, n=1, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the cn accuracy vs cna rate at a fixed number of clones (n). '''
    x = "cell_cna_rate"
    y = "cn_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.75').query('num_clones=={}'.format(n)).query('beta0==1.2')
    box_pairs = [
        ((0.05, "PERT clone"), (0.05, "PERT comp.")),
        ((0.02, "PERT clone"), (0.02, "PERT comp.")),
        ((0.00, "PERT clone"), (0.00, "PERT comp.")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across cell CNA rate\n# of clones={}'.format(n))
    ax.set_xlabel('Cell CNA rate')
    ax.set_ylabel('CN state accuracy')
    return ax


def plot_clone_effect_rep_acc(df, ax, rate=0.02, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the replication accuracy vs number of clones at a fixed cell cna rate. '''
    x = "num_clones"
    y = "rep_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.75').query('cell_cna_rate=={}'.format(rate)).query('num_clones<4').query('beta0==1.2')
    box_pairs = [
        ((3, "Kronos"), (3, "PERT comp.")),
        ((3, "PERT clone"), (3, "Kronos")),
        ((3, "PERT clone"), (3, "PERT comp.")),
        ((1, "Kronos"), (1, "PERT comp.")),
        ((1, "PERT clone"), (1, "Kronos")),
        ((1, "PERT clone"), (1, "PERT comp.")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across # of clones\nCell CNA rate={}'.format(rate))
    ax.set_xlabel('Number of clones')
    ax.set_ylabel('Replication state accuracy')
    return ax


def plot_clone_effect_cn_acc(df, ax, rate=0.02, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the CN accuracy vs number of clones at a fixed cell cna rate. '''
    x = "num_clones"
    y = "cn_accuracy"
    hue = "method"
    temp_df = df.query('alpha==10.0').query('lamb==0.75').query('cell_cna_rate=={}'.format(rate)).query('num_clones<4').query('beta0==1.2')
    box_pairs = [
        ((3, "PERT clone"), (3, "PERT comp.")),
        ((1, "PERT clone"), (1, "PERT comp.")),
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across # of clones\nCell CNA rate={}'.format(rate))
    ax.set_xlabel('Number of clones')
    ax.set_ylabel('CN state accuracy')
    return ax


def plot_alpha_effect_cn_acc(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the cn accuracy vs alpha where the hue is cell cna rate. '''
    x = "alpha"
    y = "cn_accuracy"
    hue = "method"
    temp_df = df.query('lamb==0.75').query('beta0==1.2').query('num_clones<4').query('cell_cna_rate==0.02')
    box_pairs = [
        ((10.0, 'PERT comp.'), (10.0, 'PERT clone')),
        ((5.0, 'PERT comp.'), (5.0, 'PERT clone')),
        ((15.0, 'PERT comp.'), (15.0, 'PERT clone'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across alpha')
    ax.set_ylabel('CN state accuracy')
    ax.set_xlabel('alpha')


def plot_alpha_effect_rep_acc(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the rep accuracy vs alpha where the hue is cell cna rate. '''
    x = "alpha"
    y = "rep_accuracy"
    hue = "method"
    temp_df = df.query('lamb==0.75').query('beta0==1.2').query('num_clones<4').query('cell_cna_rate==0.02')
    box_pairs = [
        ((10.0, 'PERT comp.'), (10.0, 'PERT clone')),
        ((5.0, 'PERT comp.'), (5.0, 'PERT clone')),
        ((15.0, 'PERT comp.'), (15.0, 'PERT clone')),
        ((10.0, 'PERT comp.'), (10.0, 'Kronos')),
        ((5.0, 'PERT comp.'), (5.0, 'Kronos')),
        ((15.0, 'PERT comp.'), (15.0, 'Kronos')),
        ((10.0, 'PERT clone'), (10.0, 'Kronos')),
        ((5.0, 'PERT clone'), (5.0, 'Kronos')),
        ((15.0, 'PERT clone'), (15.0, 'Kronos'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across alpha')
    ax.set_ylabel('Replication state accuracy')
    ax.set_xlabel('alpha')


def plot_lambda_effect_cn_acc(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the copy number accuracy vs lambda where the hue is cell cna rate. '''
    x = "lamb"
    y = "cn_accuracy"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0').query('num_clones==1').query('alpha==10.0').query('beta0==1.2')
    box_pairs = [
        ((0.5, 'PERT clone'), (0.5, 'PERT comp.')),
        ((0.6, 'PERT clone'), (0.6, 'PERT comp.')),
        ((0.75, 'PERT clone'), (0.75, 'PERT comp.')),
        ((0.9, 'PERT clone'), (0.9, 'PERT comp.')),
        ((0.99, 'PERT clone'), (0.99, 'PERT comp.'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_xlabel('lambda')
    ax.set_title('Sweep across lambda')
    ax.set_ylabel('CN state accuracy')


def plot_lambda_effect_rep_acc(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the replication accuracy vs lambda where the hue is the method. '''
    x = "lamb"
    y = "rep_accuracy"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0').query('num_clones==1').query('alpha==10.0').query('beta0==1.2')
    box_pairs = [
        ((0.5, 'PERT clone'), (0.5, 'PERT comp.')),
        ((0.6, 'PERT clone'), (0.6, 'PERT comp.')),
        ((0.75, 'PERT clone'), (0.75, 'PERT comp.')),
        ((0.9, 'PERT clone'), (0.9, 'PERT comp.')),
        ((0.99, 'PERT clone'), (0.99, 'PERT comp.')),
        ((0.5, 'PERT clone'), (0.5, 'Kronos')),
        ((0.6, 'PERT clone'), (0.6, 'Kronos')),
        ((0.75, 'PERT clone'), (0.75, 'Kronos')),
        ((0.9, 'PERT clone'), (0.9, 'Kronos')),
        ((0.99, 'PERT clone'), (0.99, 'Kronos')),
        ((0.5, 'Kronos'), (0.5, 'PERT comp.')),
        ((0.6, 'Kronos'), (0.6, 'PERT comp.')),
        ((0.75, 'Kronos'), (0.75, 'PERT comp.')),
        ((0.9, 'Kronos'), (0.9, 'PERT comp.')),
        ((0.99, 'Kronos'), (0.99, 'PERT comp.'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_xlabel('lambda')
    ax.set_title('Sweep across lambda')
    ax.set_ylabel('Replication state accuracy')


def plot_gc_bias_effect_cn_acc(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the phase accuracy vs GC bias coefficients in datasets with all other params fixed. '''
    x = "beta0"
    y = "cn_accuracy"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0.0').query('num_clones==1').query('alpha==10.0').query('lamb==0.75')
    box_pairs = [
        ((1.2, 'PERT comp.'), (1.2, 'PERT clone')),
        ((-1.2, 'PERT comp.'), (-1.2, 'PERT clone'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across GC bias coefficients')
    ax.set_ylabel('CN state accuracy')
    ax.set_xlabel('Beta0 (GC bias slope)')


def plot_gc_bias_effect_rep_acc(df, ax, test='t-test_ind', text_format='star', loc='inside', verbose=0, hue_order=None):
    ''' Plot the replication accuracy vs GC bias coefficients in datasets with all other params fixed. '''
    x = "beta0"
    y = "rep_accuracy"
    hue = "method"
    temp_df = df.query('cell_cna_rate==0.0').query('num_clones==1').query('alpha==10.0').query('lamb==0.75')
    box_pairs = [
        ((1.2, 'PERT comp.'), (1.2, 'PERT clone')),
        ((-1.2, 'PERT comp.'), (-1.2, 'PERT clone')),
        ((1.2, 'PERT comp.'), (1.2, 'Kronos')),
        ((-1.2, 'PERT comp.'), (-1.2, 'Kronos')),
        ((1.2, 'PERT clone'), (1.2, 'Kronos')),
        ((-1.2, 'PERT clone'), (-1.2, 'Kronos'))
    ]
    violins_with_pvals(temp_df, x, y, hue, ax, box_pairs, test=test,
                       text_format=text_format, loc=loc, verbose=verbose, hue_order=hue_order)
    ax.set_title('Sweep across GC bias coefficients')
    ax.set_ylabel('Replication state accuracy')
    ax.set_xlabel('Beta0 (GC bias slope)')


def plot_param_sweep_rep_acc(df, argv):
    """ Figure showing replication accuracy parameter sweep across num_clones and cell_cna_rate """
    fig, ax = plt.subplots(2, 5, figsize=(20, 8), tight_layout=True)

    # set the hue order to be the same for all plots
    hue_order = ['PERT comp.', 'PERT clone', 'Kronos']

    # top row shows accuracy at predicting replication states
    # show effect of varying cna rate at fixed number of clones
    plot_cna_rate_rep_acc(df, ax[0, 0], n=1, hue_order=hue_order)
    plot_cna_rate_rep_acc(df, ax[0, 1], n=3, hue_order=hue_order)
    # show effect of varying number of clones at fixed cna rates
    plot_clone_effect_rep_acc(df, ax[0, 2], rate=0.0, hue_order=hue_order)
    plot_clone_effect_rep_acc(df, ax[0, 3], rate=0.02, hue_order=hue_order)
    plot_clone_effect_rep_acc(df, ax[0, 4], rate=0.05, hue_order=hue_order)

    # merge together the two supblots in the bottom left corner
    gs = ax[1, 0].get_gridspec()
    for a in ax[1, :2]:
        a.remove()
    axbig_bottom_row = fig.add_subplot(gs[1, 0:2])

    # barplots of phase accuracies for all simulated datasets
    sns.barplot(
        data=df, x='datatag', y='rep_accuracy', hue='method', hue_order=hue_order, 
        ax=axbig_bottom_row, palette=get_methods_cmap(), saturation=1
    )
    axbig_bottom_row.set_ylabel('Replication state accuracy')
    axbig_bottom_row.set_title('All simulated datasets')

    # plot the effect of varying alpha
    plot_alpha_effect_rep_acc(df, ax[1, 2], hue_order=hue_order)

    # plot the effect of varying lambda
    plot_lambda_effect_rep_acc(df, ax[1, 3], hue_order=hue_order)
    
    # plot the effect of varying GC bias
    plot_gc_bias_effect_rep_acc(df, ax[1, 4], hue_order=hue_order)

    fig.savefig(argv.plot1, bbox_inches='tight', dpi=300)


def plot_param_sweep_cn_acc(df, argv):
    """ Figure showing copy number accuracy parameter sweep across num_clones and cell_cna_rate """
    fig, ax = plt.subplots(2, 5, figsize=(20, 8), tight_layout=True)

    # remove all rows of df where method=='Kronos' since it doesn't predict CN states
    df = df.query('method!="Kronos"')

    # hue_order
    hue_order = ['PERT comp.', 'PERT clone']

    # top row shows accuracy at predicting replication states
    # show effect of varying cna rate at fixed number of clones
    plot_cna_rate_cn_acc(df, ax[0, 0], n=1, hue_order=hue_order)
    plot_cna_rate_cn_acc(df, ax[0, 1], n=3, hue_order=hue_order)
    # show effect of varying number of clones at fixed cna rates
    plot_clone_effect_cn_acc(df, ax[0, 2], rate=0.0, hue_order=hue_order)
    plot_clone_effect_cn_acc(df, ax[0, 3], rate=0.02, hue_order=hue_order)
    plot_clone_effect_cn_acc(df, ax[0, 4], rate=0.05, hue_order=hue_order)

    # merge together the two supblots in the bottom left corner
    gs = ax[1, 0].get_gridspec()
    for a in ax[1, :2]:
        a.remove()
    axbig_bottom_row = fig.add_subplot(gs[1, 0:2])

    # barplots of phase accuracies for all simulated datasets
    sns.barplot(
        data=df, x='datatag', y='cn_accuracy', hue='method', hue_order=hue_order,
        ax=axbig_bottom_row, palette=get_methods_cmap(), saturation=1)
    axbig_bottom_row.set_ylabel('CN state accuracy')
    axbig_bottom_row.set_title('All simulated datasets')

    # plot the effect of varying alpha
    plot_alpha_effect_cn_acc(df, ax[1, 2], hue_order=hue_order)

    # plot the effect of varying lambda
    plot_lambda_effect_cn_acc(df, ax[1, 3], hue_order=hue_order)
    
    # plot the effect of varying GC bias
    plot_gc_bias_effect_cn_acc(df, ax[1, 4], hue_order=hue_order)

    fig.savefig(argv.plot2, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()

    # load the data
    df = pd.read_csv(argv.input, sep='\t')

    # rename lambda column to avoid conflict with python keyword
    df.rename(columns={'lambda': 'lamb'}, inplace=True)

    # make the first figure
    plot_param_sweep_rep_acc(df, argv)

    # make the second figure
    plot_param_sweep_cn_acc(df, argv)



if __name__ == '__main__':
    main()
