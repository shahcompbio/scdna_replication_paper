from argparse import ArgumentParser
import pandas as pd
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('input_T47D', help='long-form copy number dataframe for T47D S-phase cells with scRT data')
    p.add_argument('input_GM18507', help='long-form copy number dataframe for GM18507 S-phase cells with scRT data')
    p.add_argument('rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('output_T47D', help='scRT pseudobulks for T47D')
    p.add_argument('output_GM18507', help='scRT pseudobulks for GM18507')

    return p.parse_args()


def main():
    argv = get_args()

    # load the data for each cell line
    cn_T47D = pd.read_csv(argv.input_T47D, sep='\t')
    cn_GM18507 = pd.read_csv(argv.input_GM18507, sep='\t')

    # compute RT pseudobulks for each cell line
    bulk_rt_T47D = compute_pseudobulk_rt_profiles(cn_T47D, argv.rep_col, output_col='pseudobulk', time_col='hours', clone_col='clone_id')
    bulk_rt_GM18507 = compute_pseudobulk_rt_profiles(cn_GM18507, argv.rep_col, output_col='pseudobulk', time_col='hours', clone_col='clone_id')

    # save the pseudobulk RT profiles
    bulk_rt_T47D.to_csv(argv.output_T47D, sep='\t', index=False)
    bulk_rt_GM18507.to_csv(argv.output_GM18507, sep='\t', index=False)


if __name__ == '__main__':
    main()
