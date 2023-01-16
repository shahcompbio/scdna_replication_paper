from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('input', help='long-form copy number dataframe for G1/2-phase cells with CN data')
    p.add_argument('cn_state_col', type=str, help='column name for HMMcopy states of each bin')
    p.add_argument('dataset')
    p.add_argument('output', help='Clone cn pseudobulks')

    return p.parse_args()


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    # clone level pseudobulks
    clone_cns = compute_consensus_clone_profiles(df, argv.cn_state_col, cn_state_col=argv.cn_state_col, clone_col='clone_id')

    # add clone prefix to the clone_id
    clone_cns.columns = ['clone_{}'.format(c) for c in clone_cns.columns]
    clone_cns.reset_index(inplace=True)

    # compute sample_id column if it's not already present
    if 'sample_id' not in df.columns:
        df['sample_id'] = df['cell_id'].apply(lambda x: x.split('-')[0])

    # sample level pseudobulks
    sample_cns = compute_consensus_clone_profiles(df, argv.cn_state_col, cn_state_col=argv.cn_state_col, clone_col='sample_id')

    # add sample prefix to the clone_id
    sample_cns.columns = ['sample_{}'.format(c) for c in sample_cns.columns]
    sample_cns.reset_index(inplace=True)

    # dataset level pseudobulks
    # only compute this when the dataset_id isn't the sample_id
    if argv.dataset not in df['sample_id'].unique():
        df['dataset'] = argv.dataset
        dataset_cns = compute_consensus_clone_profiles(df, argv.cn_state_col, cn_state_col=argv.cn_state_col, clone_col='dataset')

        # add sample prefix to the clone_id
        dataset_cns.columns = ['dataset_{}'.format(c) for c in dataset_cns.columns]
        dataset_cns.reset_index(inplace=True)

        # merge clone-, sample-, and dataset-level pseudobulks into one
        bulk_cns = pd.merge(pd.merge(dataset_cns, sample_cns), clone_cns)
    else:
        # merge clone- and sample-level pseudobulks into one
        bulk_cns = pd.merge(sample_cns, clone_cns)

    bulk_cns.to_csv(argv.output, sep='\t', index=False)


if __name__ == '__main__':
    main()