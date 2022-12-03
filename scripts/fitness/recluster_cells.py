import umap
import hdbscan
import logging
import sklearn.cluster
import sklearn.mixture
import scipy.spatial
import pandas as pd
import numpy as np
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles
from scdna_replication_tools.assign_s_to_clones import assign_s_to_clones
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_g_u_in', type=str, help='filtered g1/2-phase cn data for untreated sample')
    p.add_argument('cn_g_t_in', type=str, help='filtered g1/2-phase cn data for treated sample')
    p.add_argument('cn_s_u_in', type=str, help='filtered s-phase cn data for untreated sample')
    p.add_argument('cn_s_t_in', type=str, help='filtered s-phase cn data for treated sample')
    p.add_argument('num_clusters', type=int)
    p.add_argument('cn_g_u_out', type=str, help='filtered g1/2-phase cn data for untreated sample with new clone_id added')
    p.add_argument('cn_g_t_out', type=str, help='filtered g1/2-phase cn data for treated sample with new clone_id added')
    p.add_argument('cn_s_u_out', type=str, help='filtered s-phase cn data for untreated sample with new clone_id added')
    p.add_argument('cn_s_t_out', type=str, help='filtered s-phase cn data for treated sample with new clone_id added')

    return p.parse_args()


def umap_gmm_cluster(
        cn,
        min_k=2,
        max_k=100,
    ):
    """ Cluster using gaussian mixture model and bic.
    """

    X = cn.T.values
    ks = range(min_k, max_k + 1)

    logging.info(f'trying with max k={max_k}')
    
    embedding = umap.UMAP(
        n_neighbors=15,
        min_dist=0.1,
        n_components=2,
        random_state=42,
        metric='euclidean',
    ).fit_transform(cn.fillna(0).values.T)

    gmm = []
    bics = []
    for k in ks:
        logging.info(f'trying with k={k}')
        model = sklearn.mixture.GaussianMixture(n_components=k).fit(embedding)
        bic = model.bic(embedding)
        gmm.append(model)
        bics.append(bic)

    opt_k = np.array(bics).argmax()
    logging.info(f'selected k={opt_k}')

    model = gmm[opt_k]
    labels = model.predict(embedding)

    clusters = pd.DataFrame({
        'cell_id': cn.columns, 'cluster_id': labels,
        'umap1': embedding[:, 0], 'umap2': embedding[:, 1]
    })
    
    return clusters, bics


if __name__ == '__main__':
    argv = get_args()

    # load in data from both cell cycle phases and treatment groups
    cn_g_u = pd.read_csv(argv.cn_g_u_in, sep='\t')
    cn_g_t = pd.read_csv(argv.cn_g_t_in, sep='\t')
    cn_s_u = pd.read_csv(argv.cn_s_u_in, sep='\t')
    cn_s_t = pd.read_csv(argv.cn_s_t_in, sep='\t')

    # label which cells are treated vs untreated so I can split later
    cn_g_u['treated'] = False
    cn_s_u['treated'] = False
    cn_g_t['treated'] = True
    cn_s_t['treated'] = True

    # combine both treated and untreated cells into one dataframe
    cn_g = pd.concat([cn_g_t, cn_g_u], ignore_index=True)
    cn_s = pd.concat([cn_s_t, cn_s_u], ignore_index=True)

    # use matrix of reads per million as input for clustering
    rpm_mat = pd.pivot_table(cn_g, values='rpm', columns='cell_id', index=['chr', 'start'])

    # drop loci that aren't shared amongst all cells (only a few bins)
    rpm_mat = rpm_mat.dropna()

    gmm_clusters, bics = umap_gmm_cluster(rpm_mat, min_k=argv.num_clusters, max_k=argv.num_clusters)
    cn_g = pd.merge(cn_g, gmm_clusters)
    # rename the columns
    cn_g.rename(columns={'cluster_id': 'clone_id', 'clone_id': 'sitka_clone_id'}, inplace=True)
    # convert new clone_id values from ints to letters
    cn_g['clone_id'] = cn_g['clone_id'].apply(lambda x: chr(x+65))

    ## assign the S-phase cells to these new clusters
    # compute conesensus clone profiles for assign_col
    assign_col = 'copy'
    clone_profiles = compute_consensus_clone_profiles(
        cn_g, assign_col, clone_col='clone_id', cell_col='cell_id', chr_col='chr',
        start_col='start', cn_state_col='state'
    )

    # assign S-phase cells to clones based on similarity of assign_col
    cn_s = assign_s_to_clones(
        cn_s, clone_profiles, col_name=assign_col,
        clone_col='clone_id', cell_col='cell_id', chr_col='chr', start_col='start'
    )

    # separate back into treated vs untreated
    cn_g_u_out = cn_g.loc[cn_g['treated']==False]
    cn_g_t_out = cn_g.loc[cn_g['treated']==True]
    cn_s_u_out = cn_s.loc[cn_s['treated']==False]
    cn_s_t_out = cn_s.loc[cn_s['treated']==True]

    # return one cn dataframe for each cell cycle state and treatment group
    cn_g_u_out.to_csv(argv.cn_g_u_out, sep='\t', index=False)
    cn_g_t_out.to_csv(argv.cn_g_t_out, sep='\t', index=False)
    cn_s_u_out.to_csv(argv.cn_s_u_out, sep='\t', index=False)
    cn_s_t_out.to_csv(argv.cn_s_t_out, sep='\t', index=False)

