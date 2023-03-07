from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.colors import get_cell_line_cmap
from common.plot_utils import plot_cell_cn_profile2


def get_args():
    p = ArgumentParser()

    p.add_argument('rt_bulks', help='input rt pseudobulks from multiple cell lines and model versions')
    p.add_argument('rt_diff_split', help='difference in rt values between the two cell lines for the split pyro model')
    p.add_argument('rt_diff_merged', help='difference in rt values between the two cell lines for the merged pyro model')
    p.add_argument('rt_corr', help='pairwise correlation between all rt profiles')
    p.add_argument('rt_split_chr1', help='both cell line chr1 RT profiles for split model')
    p.add_argument('rt_merged_chr1', help='both cell line chr1 RT profiles for merged model')


    return p.parse_args()


def plot_chr1_by_cell_line(df, argv):
    cell_line_cmap = get_cell_line_cmap()

    # plot chr1 bulk RT for each cell line inferred with the merged model
    fig, ax = plt.subplots(1, 1, figsize=(16,4))
    plot_cell_cn_profile2(ax, df, 'rt_merged_T47D', color=cell_line_cmap['T47D'], label='merged T47D', max_cn=None, scale_data=False, lines=True, chromosome='1')
    plot_cell_cn_profile2(ax, df, 'rt_merged_GM18507', color=cell_line_cmap['GM18507'], label='merged GM18507', max_cn=None, scale_data=False, lines=True, chromosome='1')
    ax.set_ylabel('Pseudobulk RT')
    ax.set_title('Merged PERT input')
    ax.legend()
    fig.savefig(argv.rt_merged_chr1, bbox_inches='tight', dpi=300)

    # plot chr1 bulk RT for each cell line inferred with the split model
    fig, ax = plt.subplots(1, 1, figsize=(16,4))
    plot_cell_cn_profile2(ax, df, 'rt_split_T47D', color=cell_line_cmap['T47D'], label='split T47D', max_cn=None, scale_data=False, lines=True, chromosome='1')
    plot_cell_cn_profile2(ax, df, 'rt_split_GM18507', color=cell_line_cmap['GM18507'], label='split GM18507', max_cn=None, scale_data=False, lines=True, chromosome='1')
    ax.set_ylabel('Pseudobulk RT')
    ax.set_title('Split PERT input')
    ax.legend()
    fig.savefig(argv.rt_split_chr1, bbox_inches='tight', dpi=300)


def plot_rt_corr(df, argv):
    """ plot a heatmap of all the pairwise RT correlations for each cell line when the merged and split methods are used """
    # subset the dataframe to only the columns that should be plotted
    df = df[['rt_merged_T47D', 'rt_split_T47D', 'rt_merged_GM18507', 'rt_split_GM18507']]
    # rename the columns to be more readable
    df.columns = ['T47D merged', 'T47D split', 'GM18507 merged', 'GM18507 split']
    # calculate the Pearson correlation matrix
    corr = df.corr()

    # mask the upper triangle of the heatmap
    mask = np.zeros_like(corr, dtype=bool)
    mask[np.triu_indices_from(mask)] = True

    fig = plt.figure(figsize=(6, 6))

    # plot the heatmap, including the first 2 digits of the correlation values, the mask, and whitespace between the cells
    sns.heatmap(corr, square=True, linewidths=.5, cbar_kws={"shrink": .5}, mask=mask, annot=True, fmt='.2f')
    plt.title('Pseudobulk RT correlations\nMerged vs split PERT input')

    # rotate the y-tick labels to read from left to right
    plt.yticks(rotation=0)

    fig.savefig(argv.rt_corr, bbox_inches='tight', dpi=300)


def plot_rt_diff(df, argv):
    # plot the difference in RT values between the two cell lines
    # for the merged model
    fig, ax = plt.subplots(1, 1, figsize=(12,4))
    plot_cell_cn_profile2(ax, df, 'rt_diff_merged', color='#BA0021', max_cn=None, scale_data=False, lines=True)
    ax.set_ylabel('Pseudobulk RT difference\n<--GM18507 earlier | T47D earlier -->')
    ax.set_title('Merged PERT input')
    fig.savefig(argv.rt_diff_merged, bbox_inches='tight', dpi=300)

    # for the split model
    fig, ax = plt.subplots(1, 1, figsize=(12,4))
    plot_cell_cn_profile2(ax, df, 'rt_diff_split', color='#BA0021', max_cn=None, scale_data=False, lines=True)
    ax.set_ylabel('Pseudobulk RT difference\n<--GM18507 earlier | T47D earlier -->')
    ax.set_title('Split PERT input')
    fig.savefig(argv.rt_diff_split, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()
    df = pd.read_csv(argv.rt_bulks, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # create separate plots
    plot_rt_corr(df, argv)
    plot_rt_diff(df, argv)
    plot_chr1_by_cell_line(df, argv)


if __name__ == '__main__':
    main()
