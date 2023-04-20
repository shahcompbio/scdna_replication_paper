import pandas as pd
import matplotlib.pyplot as plt
from scgenome.cnplot import plot_cell_cn_profile
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_rt_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_g', help='long-form copy number dataframe for all cells thought to be in G1/2 phase by PERT')
    p.add_argument('output', help='representative cell profiles for that were thought to be in G1/2 phase by PERT but were in S phase by flow')

    return p.parse_args()


def main():
    argv = get_args()

    # load in the input dataframe
    cn_g = pd.read_csv(argv.cn_g, sep='\t')
    cn_g['chr'] = cn_g['chr'].astype('category')

    # get the replication timing colormap
    rt_cmap = get_rt_cmap()

    fig, ax = plt.subplots(4, 2, figsize=(24, 12), tight_layout=True)

    # plot a representative T47D cell that was thought to be in G1/2 phase by PERT but was in S phase by flow
    cell_id = 'SA1044-A96139A-R03-C26'
    cell_cn = cn_g.query('cell_id=="{}"'.format(cell_id))
    
    # plot rpm colored by hmmcopy state and pert replication state
    plot_cell_cn_profile(ax[0, 0], cell_cn, 'rpm', cn_field_name='state', rawy=True)
    plot_cell_cn_profile(ax[1, 0], cell_cn, 'rpm', cn_field_name='model_rep_state', cmap=rt_cmap, rawy=True)
    # same as above but zoom in on chromosome 1
    plot_cell_cn_profile(ax[2, 0], cell_cn, 'rpm', cn_field_name='state', rawy=True, chromosome='1')
    plot_cell_cn_profile(ax[3, 0], cell_cn, 'rpm', cn_field_name='model_rep_state', cmap=rt_cmap, rawy=True, chromosome='1')
    
    # get per-cell metrics to include in the title
    frac = round(cell_cn['cell_frac_rep'].values[0], 3)
    quality = round(cell_cn['quality'].values[0], 3)
    s_prob = round(cell_cn['is_s_phase_prob'].values[0], 3)
    flow_state = cell_cn['cell_cycle_state'].values[0]
    for i in range(4):
        ax[i, 0].set_ylabel('reads per million')
        ax[i, 0].set_title('{}\nfrac rep={}, quality={}, S-phase prob={}, flow={}'.format(cell_id, frac, quality, s_prob, flow_state))
    
    # plot a representative GM18507 cell that was thought to be in G1/2 phase by PERT but was in S phase by flow
    cell_id = 'SA928-A73044A-R21-C52'
    cell_cn = cn_g.query('cell_id=="{}"'.format(cell_id))

    # plot rpm colored by hmmcopy state and pert replication state
    plot_cell_cn_profile(ax[0, 1], cell_cn, 'rpm', cn_field_name='state', rawy=True)
    plot_cell_cn_profile(ax[1, 1], cell_cn, 'rpm', cn_field_name='model_rep_state', cmap=rt_cmap, rawy=True)
    # same as above but zoom in on chromosome 1
    plot_cell_cn_profile(ax[2, 1], cell_cn, 'rpm', cn_field_name='state', rawy=True, chromosome='1')
    plot_cell_cn_profile(ax[3, 1], cell_cn, 'rpm', cn_field_name='model_rep_state', cmap=rt_cmap, rawy=True, chromosome='1')

    # get per-cell metrics to include in the title
    frac = round(cell_cn['cell_frac_rep'].values[0], 3)
    quality = round(cell_cn['quality'].values[0], 3)
    s_prob = round(cell_cn['is_s_phase_prob'].values[0], 3)
    flow_state = cell_cn['cell_cycle_state'].values[0]
    for i in range(4):
        ax[i, 1].set_ylabel('reads per million')
        ax[i, 1].set_title('{}\nfrac rep={}, quality={}, S-phase prob={}, flow={}'.format(cell_id, frac, quality, s_prob, flow_state))


    # save the figure
    fig.savefig(argv.output, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
