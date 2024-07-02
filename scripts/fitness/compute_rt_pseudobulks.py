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

    # compute sample and clone-level pseudobulk RT profiles
    bulk_rt = compute_pseudobulk_rt_profiles(df, argv.rep_col, output_col='pseudobulk', time_col='hours', clone_col='clone_id')

    # subset the dataframe to only include rows with 'U' in the datasetname as they are untreated cells
    df_u = df[df.datasetname.str.contains('U')]
    df_t = df[~df.datasetname.str.contains('U')]

    # compute the pseudobulk RT profiles for the untreated and treated cells
    bulk_rt_u = compute_pseudobulk_rt_profiles(df_u, argv.rep_col, output_col='pseudobulk', time_col='hours', clone_col='clone_id')
    bulk_rt_t = compute_pseudobulk_rt_profiles(df_t, argv.rep_col, output_col='pseudobulk', time_col='hours', clone_col='clone_id')

    # only keep the columns that are not clone-specific for treated and untreated pseudobulks
    cols = ['chr', 'start', 'pseudobulk_model_rep_state', 'pseudobulk_hours']
    bulk_rt_u = bulk_rt_u[cols]
    bulk_rt_t = bulk_rt_t[cols]

    # add the _U and _T suffixes to the columns for the untreated and treated pseudobulks
    bulk_rt_u.columns = ['chr', 'start', 'pseudobulk_model_rep_state_U', 'pseudobulk_hours_U']
    bulk_rt_t.columns = ['chr', 'start', 'pseudobulk_model_rep_state_T', 'pseudobulk_hours_T']

    # merge the untreated and treated pseudobulks
    bulk_rt_rx = pd.merge(bulk_rt_u, bulk_rt_t)

    # merge the pseudobulks with the original dataframe
    bulk_rt = pd.merge(bulk_rt, bulk_rt_rx)

    bulk_rt.to_csv(argv.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
