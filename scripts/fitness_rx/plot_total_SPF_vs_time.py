import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('s_phase_rx', help='Cell counts by clone & time for S-phase Rx')
    p.add_argument('non_s_phase_rx', help='Cell counts by clone & time for G1-phase Rx')
    p.add_argument('s_phase_unrx', help='Cell counts by clone & time for S-phase UnRx')
    p.add_argument('non_s_phase_unrx', help='Cell counts by clone & time for G1-phase UnRx')
    p.add_argument('dataset')
    p.add_argument('output_pdf', help='output pdf file fraction of each clone for S and non-S cells')

    return p.parse_args()


def compute_SPF_stdev(num_s, num_g, size=500, subsample_frac=0.75):
    # use 75% subsamples to draw 500 different SPF values (allows for error calculation)
    n = subsample_frac * (num_s + num_g)
    SPF = num_s / (num_s + num_g)
    pvals = [SPF, 1-SPF]
    a = pd.DataFrame(np.random.multinomial(n=n, pvals=pvals, size=size), columns=['subsample_num_s', 'subsample_num_g'])
    a['SPF'] = a['subsample_num_s'] / (a['subsample_num_g'] + a['subsample_num_s'])
    
    return a


def main():
    argv = get_args()

    s_rx_df = pd.read_csv(argv.s_phase_rx, sep='\t')
    g_rx_df = pd.read_csv(argv.non_s_phase_rx, sep='\t')
    s_unrx_df = pd.read_csv(argv.s_phase_unrx, sep='\t')
    g_unrx_df = pd.read_csv(argv.non_s_phase_unrx, sep='\t')

    # create new columns to combine into one df
    s_rx_df['is_s_phase'] = True
    s_rx_df['Rx'] = True
    s_unrx_df['is_s_phase'] = True
    s_unrx_df['Rx'] = False
    g_rx_df['is_s_phase'] = False
    g_rx_df['Rx'] = True
    g_unrx_df['is_s_phase'] = False
    g_unrx_df['Rx'] = False
    df = pd.concat([s_rx_df, s_unrx_df, g_rx_df, g_unrx_df], ignore_index=True)

    # compute SPF values at each timepoint and Rx status
    total_SPF_df = []
    for (t, Rx), chunk in df.groupby(['timepoint', 'Rx']):
        temp_s = chunk.query('is_s_phase==True')
        temp_g = chunk.query('is_s_phase==False')
        num_cells_s = sum(temp_s['num_cells'].values)
        num_cells_g = sum(temp_g['num_cells'].values)
        temp_df = compute_SPF_stdev(num_cells_s, num_cells_g)
        temp_df['timepoint'] = t
        temp_df['Rx'] = Rx
        temp_df['num_cells_s'] = num_cells_s
        temp_df['num_cells_g'] = num_cells_g
        total_SPF_df.append(temp_df)
    total_SPF_df = pd.concat(total_SPF_df, ignore_index=True)

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    sns.lineplot(data=total_SPF_df, x='timepoint', y='SPF', ci='sd', err_style='bars', hue='Rx', ax=ax)
    ax.set_title(argv.dataset)
    fig.savefig(argv.output_pdf, bbox_inches='tight')


if __name__ == '__main__':
    main()
