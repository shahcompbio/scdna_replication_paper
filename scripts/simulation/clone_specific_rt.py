import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('P10_input', type=str, help='pseudobulk RT profiles for simulated dataset P10')
    p.add_argument('P11_input', type=str, help='pseudobulk RT profiles for simulated dataset P11')
    p.add_argument('output', type=str, help='Heatmap of pairwise correlations between true and inferred clone RT profiles')

    return p.parse_args()

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

    # drop the 'pseudobulk_' prefix from the columns
    df.columns = df.columns.str.replace('pseudobulk_', '')

    # rename the 'model_rep_state' suffix to 'inferred'
    df.columns = df.columns.str.replace('model_rep_state', 'inferred')

    # rename the 'true_rep' suffix to 'true'
    df.columns = df.columns.str.replace('true_rep', 'true')

    # compute the pairwise correlation between all of the columns except the chr and start columns
    corr = df.corr()

    # drom the 'start' row and column from the correlation matrix
    corr = corr.drop('start', axis=0)
    corr = corr.drop('start', axis=1)

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
    fig.savefig(argv.output, bbox_inches='tight', dpi=300)

if __name__ == '__main__':
    main()
