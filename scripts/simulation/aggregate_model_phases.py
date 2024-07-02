import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-o', '--output', type=str, help='dataset ids')

    return p.parse_args()


def main():
    argv = get_args()

    s_path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/simulation/{}/s_phase_cells_pyro_composite_filtered.tsv'
    g_path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/simulation/{}/g1_phase_cells_pyro_composite_filtered.tsv'
    lq_path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/simulation/{}/model_lowqual_composite_cells.tsv'
    
    df = []
    for d in argv.datasets:
        temp_df = pd.DataFrame({
            'dataset': [d], 's_path': [s_path.format(d)],
            'g_path': [g_path.format(d)], 'lq_path': [lq_path.format(d)]
        })
        df.append(temp_df)
    df = pd.concat(df, ignore_index=True)

    df.to_csv(argv.output, sep='\t', index=False)

if __name__=='__main__':
    main()
