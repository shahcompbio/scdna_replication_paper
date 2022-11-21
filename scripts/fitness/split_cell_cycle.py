import umap
import hdbscan
import logging
import sklearn.cluster
import scipy.spatial
import pandas as pd
import numpy as np
from scdna_replication_tools.compute_consensus_clone_profiles import compute_consensus_clone_profiles
from scdna_replication_tools.assign_s_to_clones import assign_s_to_clones
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', type=str, help='filtered cn data for all cells with relevent metrics columns')
    p.add_argument('dataset')
    p.add_argument('s_out', type=str, help='cn data for S-phase cells')
    p.add_argument('g_out', type=str, help='cn data for G1/2-phase cells')

    return p.parse_args()


def compute_bic(kmeans, X):
    """ Computes the BIC metric for a given k means clustering
    Args:
        kmeans: a fitted kmeans clustering object
        X: data for which to calculate bic
    
    Returns:
        float: bic
    
    Reference: https://stats.stackexchange.com/questions/90769/using-bic-to-estimate-the-number-of-k-in-kmeans
    """
    centers = [kmeans.cluster_centers_]
    labels  = kmeans.labels_
    n_clusters = kmeans.n_clusters
    cluster_sizes = np.bincount(labels)
    N, d = X.shape

    # Compute variance for all clusters
    cl_var = (1.0 / (N - n_clusters) / d) * sum([sum(scipy.spatial.distance.cdist(X[np.where(labels == i)], [centers[0][i]], 
             'euclidean')**2) for i in range(n_clusters)])

    const_term = 0.5 * n_clusters * np.log(N) * (d+1)

    bic = np.sum([cluster_sizes[i] * np.log(cluster_sizes[i]) -
               cluster_sizes[i] * np.log(N) -
             ((cluster_sizes[i] * d) / 2) * np.log(2*np.pi*cl_var) -
             ((cluster_sizes[i] - 1) * d/ 2) for i in range(n_clusters)]) - const_term

    return bic


def kmeans_cluster(
        cn,
        min_k=2,
        max_k=100,
    ):
    """ Cluster using kmeans and bic.
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

    kmeans = []
    bics = []
    for k in ks:
        logging.info(f'trying with k={k}')
        model = sklearn.cluster.KMeans(n_clusters=k, init="k-means++").fit(embedding)
        bic = compute_bic(model, embedding)
        kmeans.append(model)
        bics.append(bic)

    opt_k = np.array(bics).argmax()
    logging.info(f'selected k={opt_k}')

    model = kmeans[opt_k]

    clusters = pd.DataFrame({
        'cell_id': cn.columns, 'cluster_id': model.labels_,
        'umap1': embedding[:, 0], 'umap2': embedding[:, 1]
    })
    
    return clusters, bics



def split_cell_cycle(cn):
    ## TODO: come up with scheme more complex than just in/out tree
    ## maybe the cells in the top quartile of ccc features within each library?
    # make sure column dtypes are correct
    cn['chr'] = cn['chr'].astype(str)
    cn['start'] = cn['start'].astype(int)
    cn['end'] = cn['end'].astype(int)
    cn['cell_id'] = cn['cell_id'].astype(str)

    cn_g = cn.loc[(cn['quality']>0.9) & (cn['corrected_breakpoints']<0)].reset_index(drop=True)
    cn_s = cn.loc[~cn['cell_id'].isin(cn_g['cell_id'].unique())].reset_index(drop=True)

    return cn_s, cn_g


if __name__ == '__main__':
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input, sep='\t')

    # split by cell cycle
    cn_s, cn_g = split_cell_cycle(cn)

    ## recluster the G1/2-phase cells using reads per million as input
    # assume the default max number of clusters is the number of sitka clones
    max_k = len(cn.clone_id.unique())
    # provide edge cases where I specify even fewer clusters
    if argv.dataset=='SA039U':
        max_k = 3
    # perform the clustering using reads per million as input
    rpm_mat = pd.pivot_table(cn_g, values='rpm', columns='cell_id', index=['chr', 'start'])
    new_clusters, bics = kmeans_cluster(rpm_mat, max_k=max_k)
    cn_g = pd.merge(cn_g, new_clusters)
    # rename the columns
    cn_g.rename(columns={'cluster_id': 'clone_id', 'clone_id': 'old_clone_id'}, inplace=True)
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
    

    # return one cn dataframe for each cell cycle state
    cn_s.to_csv(argv.s_out, sep='\t', index=False)
    cn_g.to_csv(argv.g_out, sep='\t', index=False)

