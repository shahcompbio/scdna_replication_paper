import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('s_input', type=str, help='simulated values for S-phase cells')
    p.add_argument('g1_input', type=str, help='simulated values for G1-phase cells')
    p.add_argument('s_reads', type=str, help='hmmcopy reads.csv for S-phase cells')
    p.add_argument('g1_reads', type=str, help='hmmcopy reads.csv for G1-phase cells')
    p.add_argument('s_metrics', type=str, help='hmmcopy metrics.csv for S-phase cells')
    p.add_argument('g1_metrics', type=str, help='hmmcopy metrics.csv for G1-phase cells')
    p.add_argument('s_output', type=str, help='output path for S-phase cells')
    p.add_argument('g1_output', type=str, help='output path for G1-phase cells')

    return p.parse_args()


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state'):
    ''' Compute the fraction of replicated bins for all cells in `cn` '''
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn[rep_state_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac
    return cn


def main():
    argv = get_args()

    s_input = pd.read_csv(argv.s_input, sep='\t')
    g1_input = pd.read_csv(argv.g1_input, sep='\t')
    s_reads = pd.read_csv(argv.s_reads)
    g1_reads = pd.read_csv(argv.g1_reads)
    s_metrics = pd.read_csv(argv.s_metrics)
    g1_metrics = pd.read_csv(argv.g1_metrics)

    # compute the true fraction of replicated bins for each cell
    s_input = compute_cell_frac(s_input, frac_rt_col='true_cell_frac_rep', rep_state_col='true_rep')
    g1_input['true_cell_frac_rep'] = 0.0

    # subset metrics to only include the columns we need
    metrics_columns = [
        'multiplier', 'autocorrelation_hmmcopy', 'mean_hmmcopy_reads_per_bin',
        'std_hmmcopy_reads_per_bin', 'total_mapped_reads_hmmcopy',
        'total_halfiness', 'scaled_halfiness', 'breakpoints', 'mean_copy', 'cell_id'
    ]
    s_metrics = s_metrics[metrics_columns]
    g1_metrics = g1_metrics[metrics_columns]

    # subset the hmmcopy reads columns to only include the columns we need
    reads_columns = [
        'cell_id', 'chr', 'start', 'end', 'reads', 'gc', 'map', 'copy', 'state'
    ]
    s_reads = s_reads[reads_columns]
    g1_reads = g1_reads[reads_columns]

    # merge the reads and metrics files together based on cell_id
    s_hmmcopy = s_reads.merge(s_metrics, on='cell_id')
    g1_hmmcopy = g1_reads.merge(g1_metrics, on='cell_id')

    # merge the simulated values with the hmmcopy values
    s_output = s_hmmcopy.merge(s_input)
    g1_output = g1_hmmcopy.merge(g1_input)

    # save the output
    s_output.to_csv(argv.s_output, index=False)
    g1_output.to_csv(argv.g1_output, index=False)


if __name__ == '__main__':
    main()
