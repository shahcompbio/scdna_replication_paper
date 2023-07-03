import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
from scdna_replication_tools.plot_utils import plot_cell_cn_profile2


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='table of rt profile posteriors')
    p.add_argument('--noX', action='store_true', help='binary indicator to include for version of the model where chrX is excluded')
    p.add_argument('output', type=str, help='plot of rt profiles sampled from the posterior distribution')

    return p.parse_args()


def main():
    argv = get_args()

    # read in the posterior rt profiles
    rt_profiles = pd.read_csv(argv.input)

    # set the chromosome column as a categorical variable
    rt_profiles['chr'] = rt_profiles['chr'].astype('str').astype('category')

    # plot the genome-wide learned parameters
    plot_cols = ['global_rt_0', 'ct_hgsoc_0', 'ct_tnbc_0', 'ct_htert_0', 'sig_fbi_0', 'sig_hrd_0', 'sig_td_0', 'wgd_rt_0']
    fig, ax = plt.subplots(4, 2, figsize=(12, 8), tight_layout=True, sharey=True)
    ax = ax.flatten()

    for i, col_name in enumerate(plot_cols):
        if col_name.startswith('global'):
            color = 'C0'
        elif col_name.startswith('ct'):
            color = 'C1'
        elif col_name.startswith('sig'):
            color = 'C2'
        elif col_name.startswith('wgd'):
            color = 'C3'
        plot_cell_cn_profile2(ax[i], rt_profiles, col_name, cn_field_name=None, max_cn=None,
                            chromosome=None, s=1, squashy=False, color=color, alpha=1,
                            lines=True, label=None, scale_data=False, rawy=True)
        title = col_name.replace('_0', '').replace('ct_', 'cell type = ').replace('sig_', 'signature = ').replace('wgd_rt', 'WGD')
        if argv.noX:
            title = title + ' (no chrX)'
        ax[i].set_title(title)

    fig.savefig(argv.output, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
