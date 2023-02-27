import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    parser.add_argument('input_s', help='Input tsv file for S-phase cells')
    parser.add_argument('input_g', help='Input tsv file for G1/2-phase cells')
    parser.add_argument('gc_map', help='Input csv file that contains GC and mappability of each bin')
    parser.add_argument('output', help='Output csv file that contains tracks to be used for Kronos input columns: cell_id,chr,start,end,gc_frequency,mappability,reads')
    return parser.parse_args()


def main():
    argv = get_args()

    # read input files
    df_s = pd.read_csv(argv.input_s, sep='\t')
    df_g = pd.read_csv(argv.input_g, sep='\t')

    # read GC and mappability
    df_gc = pd.read_csv(argv.gc_map, sep=',')

    # subset df_s and df_g to only include columns of interest
    df_s = df_s[['cell_id', 'chr', 'start', 'end', 'true_reads_raw']]
    df_g = df_g[['cell_id', 'chr', 'start', 'end', 'true_reads_raw']]

    # concatenate df_s and df_g
    df = pd.concat([df_s, df_g], ignore_index=True)

    # make sure chr column is string before merging with df_gc
    df['chr'] = df['chr'].astype(str)
    df_gc['chr'] = df_gc['chr'].astype(str)

    # subtract 1 from start and end to make it 0-based before merging with df_gc
    df_gc['start'] = df_gc['start'] - 1
    df = pd.merge(df, df_gc)

    # rename columns
    df = df.rename(columns={
        'true_reads_raw': 'reads',
        'gc': 'gc_frequency',
        'map': 'mappability'
    })

    # write output
    df.to_csv(argv.output, index=False)


if __name__ == '__main__':
    main()
