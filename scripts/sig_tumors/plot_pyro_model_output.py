import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scdna_replication_tools.plot_pert_output import plot_model_results
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('--cn_s', help='long-form dataframe of S-phase cells with pyro model results')
    p.add_argument('--cn_g', help='long-form dataframe of G1/2-phase cells with pyro model results')
    p.add_argument('--dataset', type=str)
    p.add_argument('--top_title_prefix', type=str)
    p.add_argument('--bottom_title_prefix', type=str)
    p.add_argument('--output', help='heatmaps of all S-phase cells sorted the same')

    return p.parse_args()



def main():
    argv = get_args()

    cn_s = pd.read_csv(argv.cn_s, dtype={'chr':str})
    cn_g = pd.read_csv(argv.cn_g, dtype={'chr':str})

    top_title_prefix = '{} {}'.format(argv.dataset, argv.top_title_prefix).replace('_', ' ')
    bottom_title_prefix = '{} {}'.format(argv.dataset, argv.bottom_title_prefix).replace('_', ' ')

    # show rpm, hmmcopy, inferred cn, inferred rep heatmaps for S-phase cells and G1/2-phase cells
    # where all the rows are sorted the same in all four heatmaps
    fig = plot_model_results(
        cn_s, cn_g, clone_col='clone_id', second_sort_col='model_tau',
        input_cn_col='state', output_cn_col='model_cn_state',
        output_rep_col='model_rep_state', rpm_col='rpm',
        rpm_title='Reads per million', input_cn_title='HMMcopy states',
        output_cn_title='PERT somatic CN states', rep_title='PERT replication states',
        top_title_prefix=top_title_prefix,
        bottom_title_prefix=bottom_title_prefix
    )

    # save the figure
    fig.savefig(argv.output, bbox_inches='tight', dpi=300)


if __name__=='__main__':
    main()
