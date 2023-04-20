import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import linregress
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser()

    parser.add_argument('-s', '--sample_counts', type=str, nargs='+', help='list of csv files with cell cycle counts (one file per sample)')
    parser.add_argument('-d', '--doubling_times', type=str, help='csv file containing doubling time and scRNA phase counts for each sample')
    parser.add_argument('-t', '--table', type=str, help='output csv file with cell cycle counts for each sample merged with doubling time and scRNA phase counts')
    parser.add_argument('-p', '--plot', type=str, help='output png file with cell cycle counts for each sample plotted against doubling time and scRNA phase counts')

    return parser.parse_args()


def main():
    argv = get_args()

    # read in the sample counts and doubling time data
    sample_counts = []
    for path in argv.sample_counts:
        temp_df = pd.read_csv(path)
        # extract the sample name from the path
        sample_name = path.split('/')[-2]
        # add a column with the sample name
        temp_df['cell_line'] = sample_name
        # add the sample counts to the list
        sample_counts.append(temp_df)
    
    # concatenate the sample counts into a single dataframe
    sample_counts = pd.concat(sample_counts, ignore_index=True)

    # convert 'dna_frac_g1' and 'dna_frac_s' to percentages
    sample_counts['pert_g1g2_pct'] = sample_counts['dna_frac_g1'] * 100
    sample_counts['pert_s_pct'] = sample_counts['dna_frac_s'] * 100

    # read in the doubling time data
    doubling_times = pd.read_csv(argv.doubling_times)

    # merge the sample counts and doubling time data
    df = sample_counts.merge(doubling_times, on='cell_line', how='inner')

    # create a list of y-axis columns to plot against 'pert_g1g2_pct'
    y_cols = ['doubling_time', 'rna_g0g1_pct']
    y_labels = ['Doubling time (h)', 'Andor scRNA G0/G1%']

    # create a panel of plots
    fig, ax = plt.subplots(nrows=1, ncols=len(y_cols), figsize=(5*len(y_cols), 5), tight_layout=True)
    ax = ax.flatten()

    # plot the data
    for i, y_col in enumerate(y_cols):
        # compute correlation coefficient between current x and y axes
        x = df['pert_g1g2_pct']
        y = df[y_col]
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        # fit a regression line to the data but don't plot the points
        sns.regplot(x='pert_g1g2_pct', y=y_col, data=df, ax=ax[i], scatter=False, color='black', line_kws={'alpha': 0.5, 'linestyle': '--'}, ci=None)
        # plot the points colored by cell line
        sns.scatterplot(x='pert_g1g2_pct', y=y_col, data=df, hue='cell_line', ax=ax[i])
        # adjust legend and axis labels
        ax[i].set_xlabel('PERT scWGS G1/2%')
        ax[i].set_ylabel(y_labels[i])
        ax[i].legend(title='Cell Line')
        ax[i].set_title('10x Gastric Cancer\nr={:.2f}, p={:.2e}'.format(r_value, p_value))
        # set the x-axis limits to be 5% wider than the data
        left_lim = max(df['pert_g1g2_pct'].min() - 5, 0)
        right_lim = min(df['pert_g1g2_pct'].max() + 5, 100)
        ax[i].set_xlim(left=left_lim, right=right_lim)
    
    # save the plot
    fig.savefig(argv.plot, dpi=300, bbox_inches='tight')

    # save the merged data
    df.to_csv(argv.table, index=False)


if __name__ == '__main__':
    main()
