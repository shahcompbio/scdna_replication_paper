from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_ccc_features import compute_ccc_features


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s_input', help='input long-form copy number dataframe for s-phase cells')
    p.add_argument('cn_g_input', help='input long-form copy number dataframe for g1-phase cells')
    p.add_argument('cn_s_out', help='same as cn_s_input with classifier features added')
    p.add_argument('cn_g_out', help='same as cn_s_input with classifier features added')


    return p.parse_args()


def main():
    argv = get_args()

    # load in data
    cn_s = pd.read_csv(argv.cn_s_input, sep='\t')
    cn_g = pd.read_csv(argv.cn_g_input, sep='\t')

    # create observed cn state column to use for computing breakpoints
    # this should serve as a proxy for hmmcopy state
    cn_s['observed_cn_state'] = cn_s['true_G1_state'] * (cn_s['true_rep'] + 1)
    cn_g['observed_cn_state'] = cn_g['true_G1_state']

    # compute the number of reads and breakpoints for each cell
    for cn in [cn_s, cn_g]:
        for cell_id, cell_cn in cn.groupby('cell_id'):
            # read count
            num_reads = sum(cell_cn['true_reads_raw'].values)
            cn.loc[cell_cn.index, 'num_reads'] = num_reads
            # breakpoints
            temp_diff = np.diff(cell_cn['observed_cn_state'].values)
            num_bk = sum(np.where(temp_diff!=0, 1, 0))
            cn.loc[cell_cn.index, 'breakpoints'] = num_bk

    # use common set of columns from cn_s and cn_g so they can be concatenated together
    keep_cols = ['chr', 'start', 'end', 'gc', 'clone_id', 'cell_id', 'true_reads_norm', 'num_reads', 'breakpoints']
    temp_cn_s = cn_s[keep_cols]
    temp_cn_g = cn_g[keep_cols]

    cn = pd.concat([temp_cn_s, temp_cn_g], ignore_index=True)

    # treat each library as a unique clone when computing the cell cycle classifier features
    cn_temp, cell_features = compute_ccc_features(
        cn, cell_col='cell_id', rpm_col='true_reads_norm', clone_col='clone_id', madn_col='madn',
        lrs_col='lrs', num_reads_col='num_reads', bk_col='breakpoints'
    )

    # merge the cell level features with the cn input
    # we don't want to use cn_temp as the output since its
    # shape might be different from the cn_input which has already been filtered
    cn_s_out = pd.merge(cn_s, cell_features)
    cn_g_out = pd.merge(cn_g, cell_features)

    # save output files
    cn_s_out.to_csv(argv.cn_s_out, sep='\t', index=False)
    cn_g_out.to_csv(argv.cn_g_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
