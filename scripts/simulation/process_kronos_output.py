import pandas as pd
import numpy as np
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser()

    parser.add_argument('kronos_input', type=str, help='Path to the Kronos input')
    parser.add_argument('kronos_output', type=str, help='Path to the Kronos output')
    parser.add_argument('out', type=str, help='Path to the merged Kronos output with proper column names')

    return parser.parse_args()


def main():
    argv = get_args()

    # Read in the Kronos input and output
    df_in = pd.read_csv(argv.kronos_input, sep='\t')
    df_out = pd.read_csv(argv.kronos_output, sep='\t')

    # make sure chr, start, end are of the correct dtypes
    df_in['chr'] = df_in['chr'].astype(str)
    df_in['start'] = df_in['start'].astype(int)
    df_in['end'] = df_in['end'].astype(int)
    df_out['chr'] = df_out['chr'].astype(str)
    df_out['start'] = df_out['start'].astype(int)
    df_out['end'] = df_out['end'].astype(int)

    # convert Rep from bool to int
    df_out['Rep'] = df_out['Rep'].astype(int)

    # add library_id column to Kronos input if it doesn't exist
    if 'library_id' not in df_in.columns:
        df_in['library_id'] = 'A'

    # subset the Kronos output to just the columns we want
    df_out = df_out[['chr', 'start', 'end', 'Cell', 'Rep', 'PercentageReplication']]

    # rename columns to match the Kronos input
    df_out.rename(columns={'Cell': 'cell_id', 'Rep': 'rt_state', 'PercentageReplication': 'frac_rt'}, inplace=True)

    # merge the Kronos input and output
    df = pd.merge(df_in, df_out, on=['chr', 'start', 'end', 'cell_id'])

    # write the merged Kronos output
    df.to_csv(argv.out, sep='\t', index=False)


if __name__ == '__main__':
    main()
