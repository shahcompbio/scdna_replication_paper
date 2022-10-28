import pandas as pd
import numpy as np
import umap
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('s_phase', help='cn input tsv file for s phase cells')
    p.add_argument('g_tree', help='cn input tsv file for g1/2 phase in the initiall tree')
    p.add_argument('g_recovered', help='cn input tsv file for g1/2 phase cells recovered by the model')
    p.add_argument('lowqual', help='cn input tsv file for low quality cells')
    p.add_argument('value_col')
    p.add_argument('dataset')
    p.add_argument('output_png', help='plot containinig umaps')

    return p.parse_args()


def main():
    argv = get_args()

    # load data
    cn_s = pd.read_csv(argv.s_phase, sep='\t')
    cn_g = pd.read_csv(argv.g_tree, sep='\t')
    cn_g_recovered = pd.read_csv(argv.g_recovered, sep='\t')
    cn_lowqual = pd.read_csv(argv.lowqual, sep='\t')

    # create column to denote cell cycle state or quality
    cn_s['phase'] = 'S'
    cn_g['phase'] = 'G1/2 tree'
    cn_g_recovered['phase'] = 'G1/2 recovered'
    cn_lowqual['phase'] = 'low quality'

    # concat into one dataframe
    cn_all = pd.concat([cn_s, cn_g, cn_g_recovered, cn_lowqual], ignore_index=True)

    # pivot to table of reads per million and create umap embedding of cells
    cn_mat = cn_all.pivot_table(index='cell_id', columns=['chr', 'start'], values=argv.value_col)
    embedding = umap.UMAP(random_state=42).fit_transform(cn_mat.values)

    # merege metric columns with embedding
    metric_cols = [
        'cell_id', 'cell_frac_rep', 'phase', 'clone_id', 'library_id',
        'total_mapped_reads_hmmcopy', 'breakpoints', 'is_s_phase_prob',
        'madn', 'lrs', 'corrected_madn', 'corrected_breakpoints', 'quality',
    ]
    metrics_df = cn_all[metric_cols].drop_duplicates()
    umap_df = pd.DataFrame({
        'cell_id': cn_mat.index, 'embedding_0': embedding[:, 0], 'embedding_1': embedding[:, 1]
    })
    umap_df = pd.merge(umap_df, metrics_df)

    # create and save the umap embeddings
    fig, ax = plt.subplots(1, 2, figsize=(10, 4), tight_layout=True)
    ax = ax.flatten()

    sns.scatterplot(data=umap_df, x='embedding_0', y='embedding_1', hue='clone_id', ax=ax[0])
    sns.scatterplot(data=umap_df, x='embedding_0', y='embedding_1', hue='phase', ax=ax[1])

    for i in range(2):
        ax[i].set_title('{} UMAP of read depth'.format(argv.dataset))
    
    fig.savefig(argv.output_png, bbox_inches='tight')


if __name__ == '__main__':
    main()
