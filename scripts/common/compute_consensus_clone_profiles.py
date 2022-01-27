from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells including input_col and clone_col')
    p.add_argument('input_col', help='column to be used for computing consesnsus profiles')
    p.add_argument('clone_col', help='column containing clone_id entries')
    p.add_argument('consensus_clones', help='consensus clone profiles used for assignment, matches input_col')

    return p.parse_args()


def main():
    argv = get_args()
    cn_g1 = pd.read_csv(argv.cn_g1, sep='\t')

    # temporarily remove columns that don't get used by infer_SPF in order to avoid
    # removing cells/loci that have NaN entries in some fields
    temp_cn_g1 = cn_g1[['cell_id', 'chr', 'start', 'end', 'state', argv.input_col, argv.clone_col]]

    # compute conesensus clone profiles
    clone_profiles = compute_consensus_clone_profiles(cn_g1, argv.input_col, clone_col=argv.clone_col)

    clone_profiles.to_csv(argv.consensus_clones, sep='\t')


if __name__ == '__main__':
    main()
