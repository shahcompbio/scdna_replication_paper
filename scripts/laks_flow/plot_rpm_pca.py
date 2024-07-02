import umap 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_phase_cmap, get_cell_line_cmap


def get_args():
    p = ArgumentParser()
    
    p.add_argument('input_s', type=str, help='filtered pyro model output for s-phase cells')
    p.add_argument('input_g', type=str, help='filtered pyro model output for g1/2-phase cells')
    p.add_argument('rpm_col', type=str, help='column name for reads per million')
    p.add_argument('output', type=str, help='figure containing UMAP embeddings colored by various features')

    return p.parse_args()


def main():
    argv = get_args()

    # load the data
    cn_s = pd.read_csv(argv.input_s, sep='\t')
    cn_g = pd.read_csv(argv.input_g, sep='\t')

    # denote the phase of each cell according to PERT
    cn_s['PERT_phase'] = 'S'
    cn_g['PERT_phase'] = 'G1/2'

    # compute the RPM matrix for all cells
    cn_all = pd.concat([cn_s, cn_g], ignore_index=True)
    cn_mat = cn_all.pivot_table(index='cell_id', columns=['chr', 'start'], values=argv.rpm_col)

    # compute the PCA embedding
    embedding = PCA(random_state=42).fit_transform(cn_mat.values)

    # save the embedding results to a dataframe
    pca_df = pd.DataFrame({
        'cell_id': cn_mat.index, 'embedding_0': embedding[:, 0], 'embedding_1': embedding[:, 1]
    })

    # create a dataframe with per-cell metrics using cn_all
    metric_cols = [
        'cell_id', 'cell_frac_rep', 'cell_cycle_state', 'PERT_phase',
        'total_mapped_reads_hmmcopy', 'breakpoints', 'is_s_phase_prob',
        'madn', 'lrs', 'corrected_madn', 'corrected_breakpoints', 'quality',
        'sample_id', 'library_id'
    ]
    metrics_df = cn_all[metric_cols].drop_duplicates()

    # rename cell_cycle_state to 'flow phase' and PERT_phase to 'PERT phase'
    metrics_df = metrics_df.rename(columns={'cell_cycle_state': 'flow phase', 'PERT_phase': 'PERT phase'})

    # create a new column for the cell line name using sample_id
    metrics_df['cell line'] = metrics_df['sample_id'].str.replace('SA928', 'GM18507').replace('SA1044', 'T47D')

    # merge the embedding and metric dataframes
    pca_df = pca_df.merge(metrics_df, on='cell_id')

    # plot the embedding colored by various features
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), tight_layout=True)
    for ax, col in zip(axes.ravel(), ['cell line', 'flow phase', 'PERT phase']):
        if col == 'cell line':
            palette = get_cell_line_cmap()
        else:
            palette = get_phase_cmap()
        sns.scatterplot(
            data=pca_df, x='embedding_0', y='embedding_1', hue=col, ax=ax, alpha=0.5, palette=palette
        )
        ax.set_title('Reads per million PCA')
        ax.set_xlabel('PC1')
        ax.set_ylabel('PC2')
    
    # save the figure
    fig.savefig(argv.output, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
