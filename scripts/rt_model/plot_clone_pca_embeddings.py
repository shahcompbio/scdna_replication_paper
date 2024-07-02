import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('embeddings', type=str, help='PCA embeddings with relevant per-profile features')
    p.add_argument('output', type=str, help='Scatterplots showing PCA embeddings colored by clone features')

    return p.parse_args()


def main():
    argv = get_args()

    # load the clone embeddings
    clone_embeddings = pd.read_csv(argv.embeddings)

    # plot the embeddings
    # only show the first 3 PCs
    fig, ax = plt.subplots(2, 3, figsize=(12, 8), tight_layout=True)
    ax = ax.flatten()

    pc1_label = 'PC1'
    pc2_label = 'PC2'
    pc3_label = 'PC3'

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC2', hue='type', size='num_cells_s', alpha=0.5, ax=ax[0])
    ax[0].set_xlabel(pc1_label)
    ax[0].set_ylabel(pc2_label)

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC2', hue='ploidy', size='num_cells_s', alpha=0.5, ax=ax[1])
    ax[1].set_xlabel(pc1_label)
    ax[1].set_ylabel(pc2_label)

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC2', hue='signature', size='num_cells_s', alpha=0.5, ax=ax[2])
    ax[2].set_xlabel(pc1_label)
    ax[2].set_ylabel(pc2_label)

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC3', hue='type', size='num_cells_s', alpha=0.5, ax=ax[3])
    ax[3].set_xlabel(pc1_label)
    ax[3].set_ylabel(pc3_label)

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC3', hue='ploidy', size='num_cells_s', alpha=0.5, ax=ax[4])
    ax[4].set_xlabel(pc1_label)
    ax[4].set_ylabel(pc3_label)

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC3', hue='signature', size='num_cells_s', alpha=0.5, ax=ax[5])
    ax[5].set_xlabel(pc1_label)
    ax[5].set_ylabel(pc3_label)

    for i in range(6):
        ax[i].set_title('Clone RT - Scores')
    
    fig.savefig(argv.output, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()