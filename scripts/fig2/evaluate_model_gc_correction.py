from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from matplotlib.colors import ListedColormap
from scipy.optimize import curve_fit


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells')
    p.add_argument('dataset')
    p.add_argument('sigma1', type=float, help='noise of read depth profiles')
    p.add_argument('gc_slope', type=float, help='slope of linear GC bias')
    p.add_argument('gc_int', type=float, help='intercept of linear GC bias')
    p.add_argument('A', type=float, help='steepness of inflection point when drawing RT state')
    p.add_argument('s_time_dist', help='distribution of S-phase cells captured (normal or uniform)')
    p.add_argument('output', help='GC bias curves before and after correction')

    return p.parse_args()


def get_rt_cmap():
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list)



def main():
    argv = get_args()
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g1 = pd.read_csv(argv.cn_g1, sep='\t')

    # set chr column to category
    cn_s.chr = cn_s.chr.astype(str)
    cn_s.chr = cn_s.chr.astype('category')
    cn_g1.chr = cn_g1.chr.astype(str)
    cn_g1.chr = cn_g1.chr.astype('category')

    fig, ax = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)
    ax = ax.flatten()

    title_second_line = 'dataset: {}, A: {}, sigma1: {}\ngc_int: {}, gc_slope: {}, s_time_dist: {}'.format(
        argv.dataset, argv.A, argv.sigma1, argv.gc_int, argv.gc_slope, argv.s_time_dist
    )

    # downsample so that not all points are plotted
    num_bins = int(1e4)
    cn_s = cn_s.sample(n=num_bins)
    cn_g1 = cn_g1.sample(n=num_bins)

    rt_cmap = get_rt_cmap()
    sns.scatterplot(data=cn_g1, x='gc', y='rpm', alpha=0.2, ax=ax[0], **{'s': 3})
    sns.scatterplot(data=cn_s, x='gc', y='rpm', hue='true_rt_state', alpha=0.2, ax=ax[1], **{'cmap': rt_cmap, 's': 3})
    sns.scatterplot(data=cn_s, x='gc', y='rpm_gc_norm', hue='true_rt_state', alpha=0.2, ax=ax[2], **{'cmap': rt_cmap, 's': 3})

    ax[0].set_title('G1/2-phase, before correction\n{}'.format(title_second_line))
    ax[1].set_title('S-phase, before correction\n{}'.format(title_second_line))
    ax[2].set_title('S-phase, after correction\n{}'.format(title_second_line))

    fig.savefig(argv.output)


if __name__ == '__main__':
    main()
