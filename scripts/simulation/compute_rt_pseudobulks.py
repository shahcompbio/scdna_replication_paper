from argparse import ArgumentParser
import pandas as pd
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('input', help='long-form copy number dataframe for S-phase cells with scRT data')
    p.add_argument('model_column', help='column to aggregate (rt state or continuous value)')
    p.add_argument('true_column', help='column to aggregate (rt state or continuous value)')
    p.add_argument('output', help='scRT pseudobulks')

    return p.parse_args()


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    model_bulk_rt = compute_pseudobulk_rt_profiles(df, argv.model_column, output_col='pseduobulk', time_col='hours', clone_col='clone_id')
    true_bulk_rt = compute_pseudobulk_rt_profiles(df, argv.true_column, output_col='true_pseduobulk', time_col='hours', clone_col='clone_id')
    bulk_rt = pd.merge(model_bulk_rt, true_bulk_rt)

    bulk_rt.to_csv(argv.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
