import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells')
    p.add_argument('metrics_input', type=str, help='metrics data for all cells')
    p.add_argument('s_out', type=str, help='cn data for S-phase cells')
    p.add_argument('g1_out', type=str, help='cn data for G1-phase cells')
    p.add_argument('g2_out', type=str, help='cn data for G2-phase cells')

    return p.parse_args()


def split_cell_cycle(cn, metrics):
    # make sure column dtypes are correct
    cn['chr'] = cn['chr'].astype(str)
    cn['start'] = cn['start'].astype(int)
    cn['end'] = cn['end'].astype(int)
    cn['cell_id'] = cn['cell_id'].astype(str)
    metrics['cell_id'] = metrics['cell_id'].astype(str)

    # subset to relevant metrics columns
    met = metrics[['cell_id', 'sample_id', 'library_id', 'cell_cycle_state', 'quality', 'total_mapped_reads_hmmcopy', 'breakpoints']]

    # subset metrics based on cell cycle state
    met_s = met.query('cell_cycle_state=="S"').reset_index(drop=True)
    met_g1 = met.query('cell_cycle_state=="G1"').reset_index(drop=True)
    met_g2 = met.query('cell_cycle_state=="G2"').reset_index(drop=True)

    # merge cn with cell cycle state specific metric dfs
    df_s = pd.merge(cn, met_s, on='cell_id')
    df_g1 = pd.merge(cn, met_g1, on='cell_id')
    df_g2 = pd.merge(cn, met_g2, on='cell_id')

    return df_s, df_g1, df_g2



if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input, sep='\t')
    metrics = pd.read_csv(argv.metrics_input, sep='\t')

    # merge metrics with cn and split by cell cycle
    df_s, df_g1, df_g2 = split_cell_cycle(cn, metrics)

    # return one cn dataframe for each cell cycle state
    df_s.to_csv(argv.s_out, sep='\t', index=False)
    df_g1.to_csv(argv.g1_out, sep='\t', index=False)
    df_g2.to_csv(argv.g2_out, sep='\t', index=False)

