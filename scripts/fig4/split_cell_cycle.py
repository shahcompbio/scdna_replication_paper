import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells with relevent metrics columns')
    p.add_argument('s_out', type=str, help='cn data for S-phase cells')
    p.add_argument('g1_out', type=str, help='cn data for G1-phase cells')
    p.add_argument('g2_out', type=str, help='cn data for G2-phase cells')

    return p.parse_args()


def split_cell_cycle(cn):
    # make sure column dtypes are correct
    cn['chr'] = cn['chr'].astype(str)
    cn['start'] = cn['start'].astype(int)
    cn['end'] = cn['end'].astype(int)
    cn['cell_id'] = cn['cell_id'].astype(str)
    
    # subset cn based on cell cycle state
    df_s = cn.query('cell_cycle_state=="S"').reset_index(drop=True)
    df_g1 = cn.query('cell_cycle_state=="G1"').reset_index(drop=True)
    df_g2 = cn.query('cell_cycle_state=="G2"').reset_index(drop=True)

    # remove low quality G1/2 phase cells
    df_g1 = df_g1.query('quality > 0.75').reset_index(drop=True)
    df_g2 = df_g2.query('quality > 0.75').reset_index(drop=True)

    # remove cells that might be in S-phase from G1/2 dfs
    df_g1 = df_g1.query('is_s_phase_prob < 0.5').reset_index(drop=True)
    df_g2 = df_g2.query('is_s_phase_prob < 0.5').reset_index(drop=True)
    df_g1 = df_g1.query('corrected_madn < 0 | corrected_breakpoints < 0').reset_index(drop=True)
    df_g2 = df_g2.query('corrected_madn < 0 | corrected_breakpoints < 0').reset_index(drop=True)

    return df_s, df_g1, df_g2



if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input, sep='\t')

    # split by cell cycle
    df_s, df_g1, df_g2 = split_cell_cycle(cn)

    # return one cn dataframe for each cell cycle state
    df_s.to_csv(argv.s_out, sep='\t', index=False)
    df_g1.to_csv(argv.g1_out, sep='\t', index=False)
    df_g2.to_csv(argv.g2_out, sep='\t', index=False)

