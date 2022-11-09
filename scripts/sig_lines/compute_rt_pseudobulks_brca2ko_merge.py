from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_SA1055', help='long-form copy number dataframe for SA1055 S-phase cells with scRT data')
    p.add_argument('cn_SA1056', help='long-form copy number dataframe for SA1056 S-phase cells with scRT data')
    p.add_argument('rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('cn_out', help='merged dataframe of the two input files')
    p.add_argument('rt_bulks', help='scRT pseudobulks')

    return p.parse_args()


def main():
    argv = get_args()

    df1 = pd.read_csv(argv.cn_SA1055, sep='\t')
    df2 = pd.read_csv(argv.cn_SA1056, sep='\t')

    # concatenate the two dataframes since they have the same columns
    df = pd.concat([df1, df2], ignore_index=True)

    # create new clone_id columns that combines both sample_id and clone_id status
    df['sample_id'] = df['cell_id'].apply(lambda x: x.split('-')[0])
    df['merged_clone_id'] = df['sample_id'] + df['clone_id']

    # compute pseudobulk rt profiles on the merged set of brca2-/- data
    bulk_rt = compute_pseudobulk_rt_profiles(df, argv.rep_col, output_col='pseduobulk', time_col='hours', clone_col='merged_clone_id')

    # return output files
    df.to_csv(argv.cn_out, sep='\t', index=False)
    bulk_rt.to_csv(argv.rt_bulks, sep='\t', index=False)


if __name__ == '__main__':
    main()
