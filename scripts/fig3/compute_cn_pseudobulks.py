from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('input', help='long-form copy number dataframe for G1/2-phase cells with CN data')
    p.add_argument('cn_state_col', type=str, help='column name for HMMcopy states of each bin')
    p.add_argument('output', help='Clone cn pseudobulks')

    return p.parse_args()


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    clone_cns = compute_consensus_clone_profiles(df, argv.cn_state_col, cn_state_col=argv.cn_state_col, clone_col='clone_id')

    clone_cns.to_csv(argv.output, sep='\t', index=True)


if __name__ == '__main__':
    main()
