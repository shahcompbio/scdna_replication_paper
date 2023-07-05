import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-cr', '--clone_rt', type=str, help='matrix of RT samples for each clone')
    p.add_argument('-f', '--features', type=str, help='clone feature metadata (e.g. signature, ploidy, cell type)')
    p.add_argument('-t', '--table', type=str, help='table mapping signature and type to each dataset')
    p.add_argument('--remove_x', action='store_true', default=False, help='filter out all loci with chr=="X"')
    p.add_argument('-e', '--embeddings', type=str, help='PCA embeddings with relevant per-profile features')
    p.add_argument('-l', '--loadings', type=str, help='PCA loadings with relevant per-loci features')
    p.add_argument('-ev', '--explained_variance_pdf', type=str, help='Curve showing explained variance by PCs')

    return p.parse_args()


def main():
    argv = get_args()

    # load the input data
    clone_features = pd.read_csv(argv.features)
    clone_rt = pd.read_csv(argv.clone_rt)
    table = pd.read_csv(argv.table, sep='\t')

    # Optional: filter out all loci with chr=='X'
    if argv.remove_x:
        clone_rt = clone_rt.query('chr!="X"')
    
    # set chr, start, end to index
    clone_rt = clone_rt.set_index(['chr', 'start', 'end'])

    num_components = 20
    pca = PCA(n_components=num_components)
    # X = StandardScaler().fit_transform(clone_rt.values.T)
    X = clone_rt.values.T
    clone_embeddings = pca.fit_transform(X)
    clone_embeddings = pd.DataFrame(clone_embeddings, index=clone_rt.columns, columns=['PC{}'.format(i+1) for i in range(num_components)])
    clone_embeddings.reset_index(inplace=True)
    clone_embeddings['dataset'] = clone_embeddings['index'].str.split('_').str[0]
    clone_embeddings['clone_id'] = clone_embeddings['index'].str.split('_').str[2].str.replace('clone', '')
    clone_embeddings.drop('index', axis=1, inplace=True)
    
    # merge with clone features
    clone_embeddings = pd.merge(clone_embeddings, clone_features, on=['dataset', 'clone_id'])
    clone_embeddings = pd.merge(clone_embeddings, table[['dataset', 'signature', 'type']], on='dataset')

    # build loadings dataframe
    clone_loadings = pd.DataFrame(pca.components_.T, index=clone_rt.index, columns=['PC{}'.format(i+1) for i in range(num_components)]).reset_index()
    clone_loadings['chr'] = clone_loadings['chr'].astype(str).astype('category')

    # save embeddings and loadings
    clone_embeddings.to_csv(argv.embeddings, index=False)
    clone_loadings.to_csv(argv.loadings, index=False)

    # plot a lineplot of the cumulative explained variance for each PC
    fig = plt.figure(figsize=(4, 4))
    sns.lineplot(x=list(range(1, num_components+1)), y=np.cumsum(pca.explained_variance_ratio_), marker='o', color='black')
    plt.xlabel('PC')
    plt.ylabel('Cumulative explained variance')
    fig.savefig(argv.explained_variance_pdf, dpi=300, bbox_inches='tight')


if __name__=='__main__':
    main()