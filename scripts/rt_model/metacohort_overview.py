import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-f', '--features', type=str, help='clone feature metadata (e.g. signature, ploidy, cell type)')
    p.add_argument('-t', '--table', type=str, help='table of input paths and metadata for each sample')
    p.add_argument('-o', '--output', type=str, help='output path for metacohort overview plot')

    return p.parse_args()


def plot_metacohort_overview(df):
    # create a multipanel figure that shows the distribution of the metadata for each dataset
    # each column is a unique clone & dataset combination
    # the first row shows the 'type' as an annotated colorbar
    # the second row shows the 'signature' as an annotated colorbar
    # the third row shows the 'condition' as an annotated colorbar
    # the fourth row shows 'ploidy' as a heatmap
    # the fifth row shows 'num_cells_s' as a barplot
    # the sixth row shows 'num_cells_g' as a barplot

    # create a figure with 6 rows and 1 column of size 12x4
    # have the first 4 rows be half the size as the last 2 rows
    fig, ax = plt.subplots(6, 1, figsize=(12, 4), sharex=True, tight_layout=True, gridspec_kw={'height_ratios': [0.2, 0.2, 0.2, 0.2, 1, 1]})

    # first row
    sns.heatmap(df['type'].map({'hTERT': 0, 'HGSOC': 1, 'TNBC': 2, 'OV2295': 3, 'T47D': 4, 'GM18507': 5}).to_frame().T, ax=ax[0], cmap='tab20', cbar=False)
    ax[0].set_ylabel('type', rotation=0, labelpad=30)

    # second row
    sns.heatmap(df['signature'].map({'NaN': 0, 'FBI': 1, 'HRD': 2, 'TD': 3}).to_frame().T, ax=ax[1], cmap='tab10', cbar=False)
    ax[1].set_ylabel('signature', rotation=0, labelpad=30)

    # third row
    sns.heatmap(df['condition'].map({'Line': 0, 'PDX': 1}).to_frame().T, ax=ax[2], cmap='tab10', cbar=False)
    ax[2].set_ylabel('condition', rotation=0, labelpad=30)

    # fourth row
    sns.heatmap(df['ploidy'].to_frame().T, ax=ax[3], cmap='viridis', cbar=False)
    ax[3].set_ylabel('ploidy', rotation=0, labelpad=30)

    # remove ytick labels from the top four rows
    for i in range(4):
        ax[i].set_yticklabels([])
        ax[i].set_yticks([])

    # fifth row
    # the x-axis is the index of the dataframe
    sns.barplot(x=df.index, y='num_cells_s', data=df, ax=ax[4], color='tab:blue')
    ax[4].set_ylabel('# S-phase\ncells')

    # sixth row
    sns.barplot(x=df.index, y='num_cells_g', data=df, ax=ax[5], color='tab:blue')
    ax[5].set_ylabel('# G1/2-phase\ncells')

    # remove the xtick labels from all rows
    for i in range(6):
        ax[i].set_xticklabels([])
        ax[i].set_xticks([])

    # set a title for the entire figure
    fig.suptitle('Metacohort Overview')

    return fig, ax


def main():
    argv = get_args()

    # load the input features and table
    # features are a csv, table is a tsv
    df = pd.read_csv(argv.features)
    table = pd.read_csv(argv.table, sep='\t')

    # drop all the one-hot-encoded columns from df that begin with type_ or signature_ 
    df = df.loc[:,~df.columns.str.startswith('type_')]
    df = df.loc[:,~df.columns.str.startswith('signature_')]

    # drop the cn_path and rt_path columns from table
    table = table.drop(['cn_path', 'rt_path', 'condition'], axis=1)

    # merge the two dataframes on the dataset column
    df = df.merge(table, on='dataset')

    # make the metacoherot overview plot
    fig, ax = plot_metacohort_overview(df)

    # save the figure
    fig.savefig(argv.output, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
