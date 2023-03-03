from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_args():
    p = ArgumentParser()

    p.add_argument('rt_table', help='table containing all the rt pseudobulk profiles for all permuted datasets')
    p.add_argument('cor_plot', help='plot showing the pairwise correlation between all permuted dataset bulk RT profiles')

    return p.parse_args()


def main():
    argv = get_args()

    # load the table containing all the rt pseudobulk profiles for all permuted datasets
    rt_wide = pd.read_csv(argv.rt_table, sep='\t')

    # remove columns that should not appear in the correlation matrix
    bad_columns = ['chr', 'start', 'end', 'gc', 'rt_split_T47D', 'rt_split_GM18507', 'rt_diff_split', 'rt_diff_merged']
    rt_wide = rt_wide.drop(columns=bad_columns)

    # rename the reference dataset columns with no permuted cells
    rt_wide = rt_wide.rename(columns={
        'rt_merged_T47D': 'T47D',
        'rt_merged_GM18507': 'GM18507',
    })

    # replace all the underscores with spaces in the column names
    rt_wide.columns = [col.replace('_', ' ') for col in rt_wide.columns]

    # compute the pairwise correlation between all bulk rt columns
    # sort the columns before computing the correlaiton to ensure the same order
    corr = rt_wide[sorted(rt_wide.columns)].corr()

    # plot the correlation matrix as a heatmap
    fig = plt.figure(figsize=(10, 10))

    # mask the upper triangle of the heatmap
    mask = np.zeros_like(corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    # plot the heatmap, including the mask and whitespace between the cells but not digit annotations
    sns.heatmap(corr, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=False)
    plt.title('Pseudobulk RT correlations\nPermuted datasets')

    # save the figure
    fig.savefig(argv.cor_plot, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
