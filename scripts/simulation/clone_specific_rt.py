import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.plot_utils import plot_cell_cn_profile2


def get_args():
    p = ArgumentParser()

    p.add_argument('P10_input', type=str, help='pseudobulk RT profiles for simulated dataset P10')
    p.add_argument('P11_input', type=str, help='pseudobulk RT profiles for simulated dataset P11')
    p.add_argument('output_corr', type=str, help='Heatmap of pairwise correlations between true and inferred clone RT profiles')
    p.add_argument('output_profiles', type=str, help='Plot the clone RT profiles directly')

    return p.parse_args()


def plot_clone_rt_profiles(df, argv):
    # plot clone RT profiles for each dataset
    # find the number of datasets
    inferred_cols = [c for c in df.columns if 'inferred' in c]
    true_cols = [c for c in df.columns if 'true' in c]

    # create a grid of subplots where each row is a dataset and the height is propotional to the number of datasets
    fig, ax = plt.subplots(2, 2, figsize=(32, 8), tight_layout=True)
    ax = ax.flatten()

    # plot the inferred RT profiles for each dataset
    for i, col in enumerate(inferred_cols):
        label = col.replace('_', ' ')
        # plot the inferred RT profile for the whole genome
        plot_cell_cn_profile2(
            ax[0], df, col, color='C{}'.format(i),
            max_cn=None, scale_data=False, lines=True, label=label
        )
        # plot the inferred RT profile for the chr1
        plot_cell_cn_profile2(
            ax[1], df, col, color='C{}'.format(i), chromosome='1',
            max_cn=None, scale_data=False, lines=True, label=label
        )
    
    # plot the true RT profiles for each dataset
    for i, col in enumerate(true_cols):
        label = col.replace('_', ' ')
        # plot the inferred RT profile for the whole genome
        plot_cell_cn_profile2(
            ax[2], df, col, color='C{}'.format(i),
            max_cn=None, scale_data=False, lines=True, label=label
        )
        # plot the inferred RT profile for the chr1
        plot_cell_cn_profile2(
            ax[3], df, col, color='C{}'.format(i), chromosome='1',
            max_cn=None, scale_data=False, lines=True, label=label
        )
    
    for a in [ax[0], ax[1]]:
        a.set_title('Inferred RT profiles')
        a.set_ylabel('RT profile\n<--late | early-->')
        a.legend(title='Clone ID', loc='upper left')
    for a in [ax[2], ax[3]]:
        a.set_title('True RT profiles')
        a.set_ylabel('RT profile\n<--late | early-->')
        a.legend(title='Clone ID', loc='upper left')
    
    # save the figure
    fig.savefig(argv.output_profiles, dpi=300, bbox_inches='tight')


def main():
    argv = get_args()

    # read in the pseudobulk RT profiles for each dataset
    df10 = pd.read_csv(argv.P10_input, sep='\t')
    df11 = pd.read_csv(argv.P11_input, sep='\t')

    # subset to just the acceptable columns in each dataset
    good_cols = [
        'chr', 'start', 
        'pseduobulk_cloneA_model_rep_state', 'true_pseduobulk_cloneA_true_rep',
        'pseduobulk_cloneB_model_rep_state', 'true_pseduobulk_cloneB_true_rep',
        'pseduobulk_cloneC_model_rep_state', 'true_pseduobulk_cloneC_true_rep',
        'pseduobulk_cloneD_model_rep_state', 'true_pseduobulk_cloneD_true_rep'
    ]

    df10 = df10[good_cols]
    df11 = df11[good_cols]

    # remove the 'true_' from the 'true_pseduobulk_' prefix from the true replication columns
    df10.columns = df10.columns.str.replace('true_pseduobulk', 'pseduobulk')
    df11.columns = df11.columns.str.replace('true_pseduobulk', 'pseduobulk')

    # correct the pseudobulk typo
    df10.columns = df10.columns.str.replace('pseduobulk', 'pseudobulk')
    df11.columns = df11.columns.str.replace('pseduobulk', 'pseudobulk')

    # create a dict that maps the clone names to the cell lines used in the simulation
    # for dataset 10
    clone_to_cell_line_10 = {
        'cloneA': 'gm06990_P10',
        'cloneB': 'gm12801_P10',
        'cloneC': 'gm12812_P10',
        'cloneD': 'gm12813_P10'
    }

    # for dataset 11
    clone_to_cell_line_11 = {
        'cloneA': 'bj_P11',
        'cloneB': 'mcf7_P11',
        'cloneC': 'hepg2_P11',
        'cloneD': 'gm12813_P11'
    }

    # use these two dicts to rename the columns in each dataset wherever the substring 'clone' appears
    df10 = df10.rename(columns={col: col.replace('clone', clone_to_cell_line_10[col.split('_')[1]]) for col in df10.columns if 'clone' in col})
    df11 = df11.rename(columns={col: col.replace('clone', clone_to_cell_line_11[col.split('_')[1]]) for col in df11.columns if 'clone' in col})

    # merge the two datasets on the chr and start columns
    df = pd.merge(df10, df11, on=['chr', 'start'])

    # convert chr column to category
    df['chr'] = df['chr'].astype('category')
    df['end'] = df['start'] + 500000 - 1

    # drop the 'pseudobulk_' prefix from the columns
    df.columns = df.columns.str.replace('pseudobulk_', '')

    # rename the 'model_rep_state' suffix to 'inferred'
    df.columns = df.columns.str.replace('model_rep_state', 'inferred')

    # rename the 'true_rep' suffix to 'true'
    df.columns = df.columns.str.replace('true_rep', 'true')

    # plot the clone RT profiles directly
    plot_clone_rt_profiles(df, argv)

    # compute the pairwise correlation between all of the columns except the chr and start columns
    corr = df.corr()

    # drop the 'start' and 'end' rows and columns from the correlation matrix
    corr = corr.drop('start', axis=0)
    corr = corr.drop('start', axis=1)
    corr = corr.drop('end', axis=0)
    corr = corr.drop('end', axis=1)

    # replace the underscores in the column and row names with spaces
    corr.columns = corr.columns.str.replace('_', ' ')
    corr.index = corr.index.str.replace('_', ' ')

    # plot the correlation matrix as a heatmap
    fig = plt.figure(figsize=(10, 10))

    # mask the upper triangle of the heatmap
    mask = np.zeros_like(corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    # plot the heatmap, including the first 2 digits of the correlation values
    sns.heatmap(corr, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='.2f')
    plt.title('Pseudobulk RT correlations\nSimulated datasets P10 and P11 with clone-specific RT')

    # save the figure
    fig.savefig(argv.output_corr, bbox_inches='tight', dpi=300)

if __name__ == '__main__':
    main()
