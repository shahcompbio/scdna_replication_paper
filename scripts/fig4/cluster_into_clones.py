import pandas as pd
from argparse import ArgumentParser
from scdna_replication_tools.cncluster import kmeans_cluster


def get_args():
    p = ArgumentParser()

    p.add_argument('g1_in', type=str, help='cn data for G1-phase cells')
    p.add_argument('g2_in', type=str, help='cn data for G2-phase cells')
    p.add_argument('value_col', help='column to use for constructing matrix from cn_g (ie rpm or state)')
    p.add_argument('dataset')
    p.add_argument('g_out', type=str, help='cn data for G1/2-phase cells with clone_id column added')

    return p.parse_args()


def cluster_into_clones(cn, min_k=2, max_k=100, clone_col='clone_id', value_col='state'):
    # convert to table where columns are cells and rows are loci
    mat = cn.pivot_table(columns='cell_id', index=['chr', 'start'], values=value_col)

    # perform kmeans clustering with using bic too pick K
    clusters = kmeans_cluster(mat, min_k=min_k, max_k=max_k)

    # change cluster_id column to clone_col
    clusters = clusters.rename(columns={'cluster_id': clone_col})

    # TODO: remove super low frequency clones (these are likely outliers)?

    # convert from numbers to letters
    clusters[clone_col] = clusters[clone_col].apply(lambda x: chr(x+65))

    cn = pd.merge(cn, clusters, on='cell_id')

    return cn


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn_g1 = pd.read_csv(argv.g1_in, sep='\t')
    cn_g2 = pd.read_csv(argv.g2_in, sep='\t')

    cn_g = pd.concat([cn_g1, cn_g2], ignore_index=True)

    if argv.dataset == 'T47D':
        min_k = 2
        max_k = 5
    elif argv.dataset == 'GM18507':
        min_k = 1
        max_k = 3
    elif argv.dataset == 'all':
        min_k = 2
        max_k = 8

    # split by cell cycle
    cn_g = cluster_into_clones(cn_g, min_k=min_k, max_k=max_k, value_col=argv.value_col)

    # return one cn dataframe for for G1 and G2 phase cells
    cn_g.to_csv(argv.g_out, sep='\t', index=False)

