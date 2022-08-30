from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_ccc_features import compute_ccc_features


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells')
    p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells including clone_id')
    p.add_argument('cn_g2', help='input long-form copy number dataframe for G2-phase cells including clone_id')
    p.add_argument('cn_s_out', help='output tsv that is same as cn_s input with features added')
    p.add_argument('cn_g1_out', help='output tsv that is same as cn_g1 input with features added')
    p.add_argument('cn_g2_out', help='output tsv that is same as cn_g2 input with features added')

    return p.parse_args()


def main():
    argv = get_args()

    # load input data from each cell cycle phase
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g1 = pd.read_csv(argv.cn_g1, sep='\t')
    cn_g2 = pd.read_csv(argv.cn_g2, sep='\t')

    # merge all phases into one dataframe
    cn = pd.concat([cn_s, cn_g1, cn_g2], ignore_index=True)

    # treat each library as a unique clone when computing the cell cycle classifier features
    cn = compute_ccc_features(
        cn, cell_col='cell_id', rpm_col='rpm', clone_col='library_id', madn_col='madn',
        lrs_col='lrs', num_reads_col='total_mapped_reads_hmmcopy', bk_col='breakpoints'
    )

    # split cn back up by the cell cycle states
    cn_s = cn.query('cell_cycle_state=="S"').reset_index(drop=True)
    cn_g1 = cn.query('cell_cycle_state=="G1"').reset_index(drop=True)
    cn_g2 = cn.query('cell_cycle_state=="G2"').reset_index(drop=True)

    # save output files
    cn_s.to_csv(argv.cn_s_out, sep='\t', index=False)
    cn_g1.to_csv(argv.cn_g1_out, sep='\t', index=False)
    cn_g2.to_csv(argv.cn_g2_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
