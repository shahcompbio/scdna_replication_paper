import pandas as pd
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()
    
    parser.add_argument('s_reads', type=str, help='hmmcopy reads.csv.gz file for S-phase cells')
    parser.add_argument('g_reads', type=str, help='hmmcopy reads.csv.gz file for G1-phase cells')
    parser.add_argument('s_metrics', type=str, help='hmmcopy metrics.csv.gz file for S-phase cells')
    parser.add_argument('g_metrics', type=str, help='hmmcopy metrics.csv.gz file for G1-phase cells')
    parser.add_argument('output_cn', type=str, help='reads.csv.gz for all cells')
    parser.add_argument('output_metrics', type=str, help='metrics.csv.gz for all cells')
    
    return parser.parse_args()


def main():
    argv = get_args()

    s_reads = pd.read_csv(argv.s_reads)
    g_reads = pd.read_csv(argv.g_reads)
    s_metrics = pd.read_csv(argv.s_metrics)
    g_metrics = pd.read_csv(argv.g_metrics)

    # concatenate reads
    reads = pd.concat([s_reads, g_reads], ignore_index=True)

    # add dummy library_id
    reads['library_id'] = 'A'

    reads.to_csv(argv.output_cn, index=False)

    # concatenate metrics
    metrics = pd.concat([s_metrics, g_metrics], ignore_index=True)

    # add dummy library_id
    metrics['library_id'] = 'A'

    metrics.to_csv(argv.output_metrics, index=False)


if __name__ == '__main__':
    main()
