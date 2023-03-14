import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()

    parser.add_argument('-c', '--clone_counts', type=str, nargs='+', help='list of csv files with cell cycle clone counts (one file per sample)')
    parser.add_argument('-d', '--doubling_times', type=str, help='csv file containing doubling time and scRNA phase counts for each sample')
    parser.add_argument('-t', '--table', type=str, help='output csv file with cell cycle counts for each sample merged with doubling time and scRNA phase counts')
    parser.add_argument('-p', '--plot', type=str, help='output png file with cell cycle counts for each sample plotted against doubling time and scRNA phase counts')

    return parser.parse_args()



def main():
    argv = get_args()

    # read in the sample counts and doubling time data
    clone_counts = []
    for path in argv.clone_counts:
        temp_df = pd.read_csv(path)
        # extract the sample name from the path
        sample_name = path.split('/')[-2]
        # add a column with the sample name
        temp_df['cell_line'] = sample_name
        # add the sample counts to the list
        clone_counts.append(temp_df)
    
    # concatenate the sample counts into a single dataframe
    clone_counts = pd.concat(clone_counts, ignore_index=True)

    # convert 'dna_frac_g1' and 'dna_frac_s' to percentages
    clone_counts['pert_g1g2_pct'] = clone_counts['dna_frac_g1'] * 100
    clone_counts['pert_s_pct'] = clone_counts['dna_frac_s'] * 100
    clone_counts['rna_frac_g1'] = clone_counts['rna_frac_g1'] * 100
    clone_counts['rna_frac_s'] = clone_counts['rna_frac_s'] * 100
    clone_counts['rna_frac_g2m'] = clone_counts['rna_frac_g2m'] * 100
    clone_counts['rna_frac_g'] = clone_counts['rna_frac_g1'] + clone_counts['rna_frac_g2m']

    # read in the doubling time data
    doubling_times = pd.read_csv(argv.doubling_times)

    # rename columns in doubling times to indicate they're cell line level measurements
    doubling_times.rename(columns={'rna_g0g1_pct': 'rna_g0g1_pct_cell_line', 'dna_g0g1_pct': 'dna_g0g1_pct_cell_line'}, inplace=True)

    # merge the sample counts and doubling time data
    df = clone_counts.merge(doubling_times, on='cell_line', how='inner')

    # rename hue and size columns for plotting purposes
    df.rename(columns={'cell_line': 'Cell Line', 'dna_num_cells_g1': '#cells'}, inplace=True)

    # create a list of y-axis columns to plot against 'pert_g1g2_pct'
    y_cols = ['doubling_time', 'rna_g0g1_pct_cell_line', 'dna_g0g1_pct_cell_line']
    y_labels = ['Doubling time (h)', 'Andor scRNA G1%', 'Andor scDNA G1/2%']

    # create a panel of plots
    fig, ax = plt.subplots(nrows=1, ncols=len(y_cols), figsize=(4*len(y_cols), 4), tight_layout=True)
    ax = ax.flatten()

    # plot the data
    for i, y_col in enumerate(y_cols):
        # fit a weighted linear regression line to the data where x='pert_g1g2_pct', y=y_col, and the weights are '#cells'
        X = sm.add_constant(df['pert_g1g2_pct'].values)
        Y = df[y_col].values
        weights = df['#cells'].values
        results = sm.WLS(Y, X, weights=weights).fit()
        # results.plot(ax=ax[i])
        # plot the points colored by cell line
        sns.scatterplot(x='pert_g1g2_pct', y=y_col, data=df, hue='Cell Line', size='#cells', ax=ax[i])
        # adjust legend and axis labels
        ax[i].set_xlabel('PERT scDNA G1/2%')
        ax[i].set_ylabel(y_labels[i])
        ax[i].set_title('10x Gastric Cancer Sublones')
        # plot the regression line
        left_lim = max(df['pert_g1g2_pct'].min() - 5, 0)
        right_lim = min(df['pert_g1g2_pct'].max() + 5, 100)
        x = np.arange(left_lim, right_lim, 1)
        y = results.params[0] + results.params[1] * x
        ax[i].plot(x, y, color='black', linestyle='--', alpha=0.5)
        # set the x-axis limits to be 5% wider than the data
        
        ax[i].set_xlim(left=left_lim, right=right_lim)
    
    # save the plot
    fig.savefig(argv.plot, dpi=300, bbox_inches='tight')

    # save the merged data
    df.to_csv(argv.table, index=False)


if __name__ == '__main__':
    main()
