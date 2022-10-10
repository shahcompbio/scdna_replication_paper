import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-o', '--output', type=str, help='dataset ids')

    return p.parse_args()


def main():
    argv = get_args()

    bulk_path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/fig2/{}/s_phase_cells_bulk_filtered.tsv'
    pyro_path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/fig2/{}/s_phase_cells_pyro_filtered.tsv'
    comp_path = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/analysis/fig2/{}/s_phase_cells_pyro_composite_filtered.tsv'
    
    df = []
    for d in argv.datasets:
        temp_df = pd.DataFrame({
            'dataset': [d], 'bulk_path': [bulk_path.format(d)],
            'clone_path': [pyro_path.format(d)], 'comp_path': [comp_path.format(d)]
        })
        df.append(temp_df)
    df = pd.concat(df, ignore_index=True)

    df.to_csv(argv.output, sep='\t', index=False)

if __name__=='__main__':
    main()
