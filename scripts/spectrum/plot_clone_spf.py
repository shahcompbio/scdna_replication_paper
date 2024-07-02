import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from argparse import ArgumentParser
from scdna_replication_tools.plot_utils import get_clone_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='cell cycle clone counts tsv')
    p.add_argument('dataset', type=str, help='name of this dataset')
    p.add_argument('plot', type=str, help='figure showing the distribution of cell cycle phases for each clone')

    return p.parse_args()
    

def plot_clone_spf(df, argv):
    ''' Plot the distribution of cell cycle phases for each clone '''
    fig, ax = plt.subplots(1, 3, figsize=(15, 5), tight_layout=True)

    pthresh = 1e-2

    # create custom legend for clones & libraries
    clone_cmap = get_clone_cmap()
    library_cmap = {}
    clone_legend_elements = [
        Line2D([0], [0], marker='^', color='w', label='enriched', markerfacecolor='k', markersize=10),
        Line2D([0], [0], marker='v', color='w', label='depleted', markerfacecolor='k', markersize=10)
    ]
    library_legend_elements = clone_legend_elements.copy()
    for i, c in enumerate(sorted(df.clone_id.unique())):
        color = clone_cmap[c]
        clone_legend_elements.append(Patch(facecolor=color, label=c))

    for i, l in enumerate(df.library_id.unique()):
        if l == 'all':
            continue  # skip the 'all' library
        color = 'C{}'.format(i)
        library_cmap[l] = color
        library_legend_elements.append(Patch(facecolor=color, label=l))
    
    # the left panel shows the distribution of cell cycle phases for each clone where library_id=='all'
    df2 = df.query('library_id=="all"').reset_index(drop=True)
    for i, row in df2.iterrows():
        clone_id = row['clone_id']
        if row['positive_p_adj'] < pthresh:
            ax[0].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id], marker='^')
        elif row['negative_p_adj'] < pthresh:
            ax[0].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id], marker='v')
        else:
            ax[0].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id])

    # draw scatterplot comparing the relative fraction of each clone in S vs G1/2 phases
    # each point corresponds to a unique clone+library combination
    df2 = df.query('library_id!="all"').reset_index(drop=True)
    for i, row in df2.iterrows():
        clone_id = row['clone_id']
        library_id = row['library_id']
        if row['positive_p_adj'] < pthresh:
            ax[1].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id], marker='^')
            ax[2].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=library_cmap[library_id], marker='^')
        elif row['negative_p_adj'] < pthresh:
            ax[1].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id], marker='v')
            ax[2].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=library_cmap[library_id], marker='v')
        else:
            ax[1].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=clone_cmap[clone_id])
            ax[2].scatter(x=row['clone_frac_g'], y=row['clone_frac_s'], c=library_cmap[library_id])

    # draw y=x line where we expect "neutral" clones to lie
    # limits for the left panel
    lims = [
        np.min([ax[0].get_xlim(), ax[0].get_ylim()]),  # min of both axes
        np.max([ax[0].get_xlim(), ax[0].get_ylim()]),  # max of both axes
    ]
    ax[0].plot(lims, lims, 'k--', alpha=0.25, zorder=0)
    # the limits are the same for the right two panels
    lims = [
        np.min([ax[1].get_xlim(), ax[1].get_ylim()]),  # min of both axes
        np.max([ax[1].get_xlim(), ax[1].get_ylim()]),  # max of both axes
    ]
    ax[1].plot(lims, lims, 'k--', alpha=0.25, zorder=0)
    ax[2].plot(lims, lims, 'k--', alpha=0.25, zorder=0)

    for i in range(2):
        ax[i].legend(handles=clone_legend_elements, title='Clone ID')
        ax[i].set_xlabel('Fraction of G1/2 cells in library assigned to clone')
        ax[i].set_ylabel('Fraction of S cells in library assigned to clone')
        ax[i].set_title('{}\nRelative proliferation rate of clones'.format(argv.dataset))
    
    ax[2].legend(handles=library_legend_elements, title='Library ID')
    ax[2].set_xlabel('Fraction of G1/2 cells in library assigned to clone')
    ax[2].set_ylabel('Fraction of S cells in library assigned to clone')
    ax[2].set_title('{}\nRelative proliferation rate of clones'.format(argv.dataset))

    # save spf figure
    fig.savefig(argv.plot, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()
    
    # load table of cell cycle clone counts
    df = pd.read_csv(argv.input)

    # plot S-phase fractions at the clone and sample levels
    # save both the plot and table used to make the plot
    plot_clone_spf(df, argv)


if __name__ == '__main__':
    main()

