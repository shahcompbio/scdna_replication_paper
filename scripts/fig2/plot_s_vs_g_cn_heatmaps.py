import matplotlib
matplotlib.use('Agg')
import sys
import scgenome.cnplot
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from scgenome import cncluster
from matplotlib.patches import Patch


def get_args():
    p = ArgumentParser()

    p.add_argument('s_phase', help='input cell_id to clone_id tsv file for s phase cells')
    p.add_argument('non_s_phase', help='input cell_id to clone_id tsv file for non s phase cells')
    p.add_argument('value_col')
    p.add_argument('dataset')
    p.add_argument('output_pdf', help='output pdf file fraction of each clone for S and non-S cells')

    return p.parse_args()


def plot_cn_heatmap(cn_g, cn_s, figsize=(18,9), dataset=None, value_col='state'):
    ''' Plot clustered cell cn matrices for S- and G1-phase copy number states with no clones labeled '''
    fig, ax = plt.subplots(1, 2, figsize=figsize, tight_layout=True)
    ax = ax.flatten()

    plot_data_g1 = scgenome.cnplot.plot_clustered_cell_cn_matrix(
        ax[0], cn_g, value_col
    )

    plot_data_s = scgenome.cnplot.plot_clustered_cell_cn_matrix(
        ax[1], cn_s, value_col
    )
    
    if dataset:
        ax[0].set_title('{}: G1/2-phase'.format(dataset))
        ax[1].set_title('{}: S-phase'.format(dataset))
    

    return fig


def main():
    argv = get_args()

    cn_s = pd.read_csv(argv.s_phase, sep='\t')
    cn_g = pd.read_csv(argv.non_s_phase, sep='\t')

    cn_s.chr = cn_s.chr.astype(str)
    cn_g.chr = cn_g.chr.astype(str)

    cn_g['cluster_id'] = 'A'
    cn_s['cluster_id'] = 'A'

    fig = plot_cn_heatmap(cn_g, cn_s, dataset=argv.dataset, value_col=argv.value_col)

    fig.savefig(argv.output_pdf, bbox_inches='tight')


if __name__ == '__main__':
    main()
