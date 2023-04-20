import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
from sklearn.decomposition import PCA
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_phase_cmap, get_clone_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('s_phase', help='cn input tsv file for s phase cells')
    p.add_argument('g_phase', help='cn input tsv file for g1/2 phase cells')
    p.add_argument('lowqual', help='cn input tsv file for low quality cells')
    p.add_argument('value_col')
    p.add_argument('dataset')
    p.add_argument('output_png', help='plot containinig pca embeddings')

    return p.parse_args()


def main():
    argv = get_args()

    # load data
    cn_s = pd.read_csv(argv.s_phase, sep='\t')
    cn_g = pd.read_csv(argv.g_phase, sep='\t')
    cn_lowqual = pd.read_csv(argv.lowqual, sep='\t')

    # create column to denote cell cycle state or quality
    cn_s['PERT_phase'] = 'S'
    cn_g['PERT_phase'] = 'G1/2'
    cn_lowqual['PERT_phase'] = 'LQ'

    # concat into one dataframe
    cn_all = pd.concat([cn_s, cn_g, cn_lowqual], ignore_index=True)

    # pivot to table of reads per million and create umap embedding of cells
    cn_mat = cn_all.pivot_table(index='cell_id', columns=['chr', 'start'], values=argv.value_col)
    embedding = PCA(random_state=42).fit_transform(cn_mat.values)

    # merege metric columns with embedding
    metric_cols = [
        'cell_id', 'cell_frac_rep', 'PERT_phase', 'clone_id', 'library_id',
        'total_mapped_reads_hmmcopy', 'breakpoints', 'is_s_phase_prob',
        'madn', 'lrs', 'corrected_madn', 'corrected_breakpoints', 'quality',
    ]
    metrics_df = cn_all[metric_cols].drop_duplicates()
    pca_df = pd.DataFrame({
        'cell_id': cn_mat.index, 'embedding_0': embedding[:, 0], 'embedding_1': embedding[:, 1]
    })
    pca_df = pd.merge(pca_df, metrics_df)

    # rename clone_id to 'Clone ID' and 'PERT_phase' to 'PERT phase' for plotting
    pca_df.rename(columns={'clone_id': 'Clone ID', 'PERT_phase': 'PERT phase'}, inplace=True)

    # create and save the pca embeddings
    fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True)
    ax = ax.flatten()

    phase_cmap = get_phase_cmap()
    clone_cmap = get_clone_cmap()
    clone_order = sorted(pca_df['Clone ID'].unique())

    sns.scatterplot(data=pca_df, x='embedding_0', y='embedding_1', hue='Clone ID', alpha=0.5, ax=ax[0], palette=clone_cmap, hue_order=clone_order)
    sns.scatterplot(data=pca_df, x='embedding_0', y='embedding_1', hue='PERT phase', alpha=0.5, ax=ax[1], palette=phase_cmap)

    for i in range(2):
        ax[i].set_title('{} reads per million PCA'.format(argv.dataset))
        ax[i].set_xlabel('PC1')
        ax[i].set_ylabel('PC2')
    
    fig.savefig(argv.output_png, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
