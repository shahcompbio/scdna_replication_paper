from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_ccc_features import compute_ccc_features


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s_input', help='input long-form copy number dataframe for s-phase cells')
    p.add_argument('cn_g_input', help='input long-form copy number dataframe for g1-phase cells')
    p.add_argument('laks_predictions', help='predicted cell cycle phase from Laks et al classifier')
    p.add_argument('cn_s_out', help='same as cn_s_input with classifier features added')
    p.add_argument('cn_g_out', help='same as cn_s_input with classifier features added')


    return p.parse_args()


def main():
    argv = get_args()

    # load in data
    cn_s = pd.read_csv(argv.cn_s_input)
    cn_g = pd.read_csv(argv.cn_g_input)
    laks_phases = pd.read_csv(argv.laks_predictions)

    # use common set of columns from cn_s and cn_g so they can be concatenated together
    keep_cols = ['chr', 'start', 'end', 'gc', 'clone_id', 'cell_id', 'reads', 'total_mapped_reads_hmmcopy', 'breakpoints']
    temp_cn_s = cn_s[keep_cols]
    temp_cn_g = cn_g[keep_cols]

    cn = pd.concat([temp_cn_s, temp_cn_g], ignore_index=True)

    # treat each library as a unique clone when computing the cell cycle classifier features
    cn_temp, cell_features = compute_ccc_features(
        cn, cell_col='cell_id', rpm_col='reads', clone_col='clone_id', madn_col='madn',
        lrs_col='lrs', num_reads_col='total_mapped_reads_hmmcopy', bk_col='breakpoints'
    )

    # add column to laks_phases to indicate if the cell is in s-phase
    # if is_s_phase is True, then the cell is in 'S' phase, otherwise it is in 'G1/2' phase
    laks_phases['laks_phase'] = laks_phases['is_s_phase'].apply(lambda x: 'S' if x else 'G1/2')
    # add laks_ prefix to is_s_phase and is_s_phase_prob columns in laks_phases
    laks_phases = laks_phases.rename(columns={'is_s_phase': 'laks_is_s_phase', 'is_s_phase_prob': 'laks_is_s_phase_prob'})
    # merge the laks predictions with cell_features
    cell_features = pd.merge(cell_features, laks_phases)
    
    # merge the cell level features with the cn input
    # we don't want to use cn_temp as the output since its
    # shape might be different from the cn_input which has already been filtered
    cn_s_out = pd.merge(cn_s, cell_features)
    cn_g_out = pd.merge(cn_g, cell_features)

    # save output files
    cn_s_out.to_csv(argv.cn_s_out, index=False)
    cn_g_out.to_csv(argv.cn_g_out, index=False)


if __name__ == '__main__':
    main()
