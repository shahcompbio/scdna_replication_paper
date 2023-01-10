from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='long-form scRT results every dataset')
    p.add_argument('--table', help='table of cell time in S-phase results across the whole cohort')
    p.add_argument('--plot', help='histogram of cell time in S-phase across the whole cohort')

    return p.parse_args()


def load_data(argv):
    """ Read in the scRT for each data cohort, only keeping the relevant per-cell columns. """
    df = []

    # loop through the input files
    for path in argv.input:
        # read in the tsv
        temp_df = pd.read_csv(path, sep='\t')
        # use the path name to extract the dataset name
        d = path.split('/')[2]
        # add the dataset name as a column
        temp_df['dataset'] = d
        # keep only the relevant columns
        temp_df = temp_df[['dataset', 'cell_id', 'library_id', 'clone_id', 'cell_frac_rep', 'model_tau']].drop_duplicates().reset_index(drop=True)
        # append to the list
        df.append(temp_df)
    
    # concatenate the list of dataframes into one
    df = pd.concat(df, ignore_index=True)

    return df


def main():
    argv = get_args()

    # load the data across the entire cohort
    df = load_data(argv)

    fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True)
    ax = ax.flatten()

    # violinplot
    sns.histplot(data=df, x='cell_frac_rep', hue='dataset', ax=ax[0])
    sns.kdeplot(data=df, x='cell_frac_rep', hue='dataset', ax=ax[1])

    for i in range(2):
        ax[i].set_xlabel('Inferred fraction of replicated bins')
        ax[i].set_title('Distribution of cells\nwithin S-phase')
    
    fig.savefig(argv.plot, dpi=300, bbox_inches='tight')



if __name__ == '__main__':
    main()