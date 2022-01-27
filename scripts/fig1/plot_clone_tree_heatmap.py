import os
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dst
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('clone_states', help='input clone copynumber states tsv file')
    p.add_argument('dataset')
    p.add_argument('output_pdf', help='output pdf file of dendrogram next to heatmap of consensus clone profiles')

    return p.parse_args()


def main():
    argv = get_args()
    cn = pd.read_csv(argv.clone_states, sep='\t')

    cn.set_index(['chr', 'start', 'end'], inplace=True)

    # compute linkage (input for dendrogram)
    D = dst.squareform(dst.pdist(cn.T, 'cityblock'))
    Y = sch.linkage(D, method='complete')

    # convert cn to long-form df (input for heatmap)
    cn2 = cn.reset_index().melt(id_vars=['chr', 'start', 'end'], var_name='cell_id', value_name='state')
    cn2['cluster_id'] = 0

    # create figure with heatmap being 10x wider than dendrogram
    fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [1, 10]}, figsize=(11, 4), tight_layout=True)
    ax = ax.flatten()

    # plot dendrogram
    Z = sch.dendrogram(Y, ax=ax[0], color_threshold=-1, labels=list(cn.columns), orientation='left')
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['left'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['bottom'].set_visible(False)
    ax[0].tick_params(top=False, bottom=False, left=False, right=False, labelright=True, labelbottom=False)

    # use dendrogram leaves to set the order each clone should appear in the heatmap
    cn2['dendrogram_order'] = 0
    i = 0
    for clone_id, chunk in cn2.groupby('cell_id'):
        cn2.loc[chunk.index, 'dendrogram_order'] = np.where(np.array(Z['leaves'])==i)[0][0]
        i += 1

    # plot the consensus clone heatmap
    plot_clustered_cell_cn_matrix(ax[1], cn2, 'state', secondary_field_name='dendrogram_order')
    ax[1].tick_params(labelleft=False)
    ax[1].set_title('{} clone pseudobulks'.format(argv.dataset))

    fig.savefig(argv.output_pdf, bbox_inches='tight')


if __name__ == '__main__':
    main()
