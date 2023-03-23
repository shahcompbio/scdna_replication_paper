from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_htert_cmap


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='Twidth value tables for every dataset')
    p.add_argument('-d', '--dataset', type=str, nargs='+')
    p.add_argument('-l', '--cell_type_labels', type=str, nargs='+')
    p.add_argument('--table', help='table of all the computed t_width values')
    p.add_argument('--plot', help='Summary plots of Tw values across all datasets')

    return p.parse_args()


def load_data(argv):
    # load long-form cn data for each dataset passed into the model
    i = 0
    df = []
    for path in argv.input:
        temp_df = pd.read_csv(path, sep='\t')
        temp_df['dataset'] = argv.dataset[i]
        temp_df['cell_type'] = argv.cell_type_labels[i]
        df.append(temp_df)
        i += 1

    # concatenate into one df
    df = pd.concat(df, ignore_index=True)
    return df


def main():
    argv = get_args()
    
    # load all the s-phase cells from all datasets into one dataframe
    df = load_data(argv)

    print(df)

    # save figure of twidth curves
    fig, ax = plt.subplots(2, 2, figsize=(8, 8), tight_layout=True)
    ax = ax.flatten()

    htert_cmap = get_htert_cmap()

    # barplots of T-width values
    # 1. T-width values for each genetic condition
    sns.barplot(data=df.query('per_cell==False'), x='cell_type', y='T-width', ax=ax[0], palette=htert_cmap)
    ax[0].tick_params(axis='x', labelrotation=45)
    ax[0].set_title('Cellular scRT heterogeneity')
    ax[0].set_xlabel('Genetic background (cell line)')

    # 2. T-width values for each dataset
    sns.barplot(data=df.query('per_cell==False'), x='dataset', y='T-width', ax=ax[1], palette=htert_cmap)
    ax[1].tick_params(axis='x', labelrotation=45)
    ax[1].set_title('Cellular scRT heterogeneity')
    ax[1].set_xlabel('')

    # compute pearson correlation coefficient between T-width and number of clones
    x = df.query('per_cell==False')['num_clones']
    y = df.query('per_cell==False')['T-width']
    slope, intercept, r_value, p_value, std_err = linregress(x, y)

    # scatterplots of T-width values vs number of cells
    # 3. T-width values for each genetic condition
    sns.regplot(data=df.query('per_cell==False'), x='num_clones', y='T-width', ax=ax[2], scatter=False, color='black', line_kws={'linestyle':'--'})
    sns.scatterplot(data=df.query('per_cell==False'), x='num_clones', y='T-width', hue='cell_type', ax=ax[2], palette=htert_cmap)
    ax[2].set_title('Cellular scRT heterogeneity\nPearson r={:.2f}, p={:.2e}'.format(r_value, p_value))
    ax[2].set_xlabel('Number of clones')
    # change xlim to be -1 from the current left limit and +1 from the current right limit
    ax[2].set_xlim(ax[2].get_xlim()[0]-1, ax[2].get_xlim()[1]+1)

    # 4. T-width values for each dataset
    # compute linear regression line between T-width and number of clones
    sns.regplot(data=df.query('per_cell==False'), x='num_clones', y='T-width', ax=ax[3], scatter=False, color='black', line_kws={'linestyle':'--'})
    sns.scatterplot(data=df.query('per_cell==False'), x='num_clones', y='T-width', hue='dataset', ax=ax[3], palette=htert_cmap)    
    ax[3].set_title('Cellular scRT heterogeneity\nPearson r={:.2f}, p={:.2e}'.format(r_value, p_value))
    ax[3].set_xlabel('Number of clones')
    # change xlim to be -1 from the current left limit and +1 from the current right limit
    ax[3].set_xlim(ax[3].get_xlim()[0]-1, ax[3].get_xlim()[1]+1)

    
    fig.savefig(argv.plot, bbox_inches='tight', dpi=300)

    # save a table of all the computed T-width values
    df.to_csv(argv.table, sep='\t', index=False)


if __name__ == '__main__':
    main()
