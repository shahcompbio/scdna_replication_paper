from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('input', help='long-form copy number dataframe for S-phase cells with scRT data')
    p.add_argument('rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('output', help='scRT pseudobulks')

    return p.parse_args()


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    bulk_rt = compute_pseudobulk_rt_profiles(df, argv.rep_col, output_col='pseduobulk', time_col='hours', clone_col='clone_id')

    bulk_rt.to_csv(argv.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
