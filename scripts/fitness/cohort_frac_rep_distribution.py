from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='long-form scRT results every dataset')
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
        temp_df = temp_df[['dataset', 'label', 'datasetname', 'cell_id', 'library_id', 'clone_id', 'cell_frac_rep', 'model_tau']].drop_duplicates().reset_index(drop=True)
        # add a 'cisplatin' column that is False if a 'U' appears in the datasetname of the cell, True otherwise
        temp_df['cisplatin'] = temp_df['datasetname'].apply(lambda x: 'Rx-' if 'U' in x else 'Rx+')
        
        # append to the list
        df.append(temp_df)
    
    # concatenate the list of dataframes into one
    df = pd.concat(df, ignore_index=True)

    return df


def main():
    argv = get_args()

    # load the data across the entire cohort
    df = load_data(argv)

    dataset_list = df['dataset'].unique()

    nrows = 2
    ncols = 3 + len(dataset_list)
    fig, ax = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows), tight_layout=True)

    # plot kdeplots of the cell S-phase times
    # top row is with common_norm=True, bottom row is with common_norm=False
    # plot a curve for every experiment (i.e. label)
    sns.kdeplot(data=df, x='cell_frac_rep', hue='label', ax=ax[0, 0])
    sns.kdeplot(data=df, x='cell_frac_rep', hue='label', common_norm=False, ax=ax[1, 0])

    # plot curves for on- vs off-cisplatin
    sns.kdeplot(data=df, x='cell_frac_rep', hue='cisplatin', ax=ax[0, 1])
    sns.kdeplot(data=df, x='cell_frac_rep', hue='cisplatin', common_norm=False, ax=ax[1, 1])

    # plot curves for each dataset
    sns.kdeplot(data=df, x='cell_frac_rep', hue='dataset', ax=ax[0, 2])
    sns.kdeplot(data=df, x='cell_frac_rep', hue='dataset', common_norm=False, ax=ax[1, 2])

    # set the x-axis labels and titles
    for i in range(nrows):
        for j in range(3):
            ax[i, j].set_xlabel('Inferred fraction of replicated loci')
            ax[i, j].set_title('Cell S-phase times')

    # plot on vs off cisplatin curves for each dataset
    j = 3
    for d in dataset_list:
        sns.kdeplot(data=df.query('dataset=="{}"'.format(d)), x='cell_frac_rep', hue='cisplatin', ax=ax[0, j])
        sns.kdeplot(data=df.query('dataset=="{}"'.format(d)), x='cell_frac_rep', hue='cisplatin', common_norm=False, ax=ax[1, j])
        for i in range(nrows):
            ax[i, j].set_xlabel('Inferred fraction of replicated loci')
            ax[i, j].set_title('{} cell S-phase times'.format(d))
        j += 1
    
    
    fig.savefig(argv.plot, dpi=300, bbox_inches='tight')



if __name__ == '__main__':
    main()