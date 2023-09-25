import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['pdf.use14corefonts'] = True
import seaborn as sns
from scdna_replication_tools.plot_utils import plot_cell_cn_profile2
from matplotlib.colors import LinearSegmentedColormap
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('features', type=str, help='csv table of clone features')
    p.add_argument('metacohort', type=str, help='tsv table of metacohort annotations')
    p.add_argument('pca_embeddings', type=str, help='csv table of PCA embeddings')
    p.add_argument('pca_regression', type=str, help='csv table of PCA regression results')
    p.add_argument('rt_profiles', type=str, help='csv table of rt profile posteriors')
    p.add_argument('output', type=str, help='plot that should be used for main figure')

    return p.parse_args()


cell_type_cdict = {
    'hTERT': 'lightsteelblue', 0: 'lightsteelblue',
    'HGSOC': 'teal', 1: 'teal',
    'TNBC': 'salmon', 2: 'salmon',
    'OV2295': 'lightgreen', 3: 'lightgreen',
    'T47D': 'orchid', 4: 'orchid',
    'GM18507': 'khaki', 5: 'khaki',
}
cell_type_cmap = LinearSegmentedColormap.from_list('cell_type_cmap', list(cell_type_cdict.values()), N=len(cell_type_cdict))

signature_cdict = {
    'FBI': 'plum', 0: 'plum',
    'HRD': 'cyan', 1: 'cyan',
    'TD': 'coral', 2: 'coral',
}
signature_cmap = LinearSegmentedColormap.from_list('signature_cmap', list(signature_cdict.values()), N=len(signature_cdict))

condition_cdict = {
    'Line': 'tan', 0: 'tan',
    'PDX': 'lightskyblue', 1: 'lightskyblue',
}
condition_cmap = LinearSegmentedColormap.from_list('condition_cmap', list(condition_cdict.values()), N=len(condition_cdict))

ploidy_cdict = {2:'#CCCCCC', 3:'#FDCC8A', 4:'#FC8D59', 5:'#E34A33'}
ploidy_cmap = LinearSegmentedColormap.from_list('ploidy_cmap', list(ploidy_cdict.values()), N=len(ploidy_cdict))

feature_cdict = {'ploidy': 'C3', 'type': 'C1', 'signature': 'C2'}


def plot_metacohort_overview(ax, df):
    # first row
    # plot the type as a colorbar where each value corresponds to a different color
    ax[0,0].imshow(df['type'].map({'hTERT': 0, 'HGSOC': 1, 'TNBC': 2, 'OV2295': 3, 'T47D': 4, 'GM18507': 5}).to_frame().T, cmap=cell_type_cmap) 
    # set a y-axis label
    ax[0,0].set_ylabel('type', rotation=0, labelpad=30)

    # second row
    ax[1,0].imshow(df['signature'].map({'NaN': 0, 'FBI': 1, 'HRD': 2, 'TD': 3}).to_frame().T, cmap=signature_cmap)
    ax[1,0].set_ylabel('signature', rotation=0, labelpad=30)

    # third row
    ax[2,0].imshow(df['condition'].map({'Line': 0, 'PDX': 1}).to_frame().T, cmap=condition_cmap)
    ax[2,0].set_ylabel('condition', rotation=0, labelpad=30)

    # fourth row
    ax[3,0].imshow(df['ploidy'].astype(int).to_frame().T, cmap=ploidy_cmap)
    ax[3,0].set_ylabel('ploidy', rotation=0, labelpad=30)

    # fifth row
    # the x-axis is the index of the dataframe
    sns.barplot(x=df.index, y='num_cells_s', data=df, ax=ax[4,0], color='tab:blue')
    ax[4,0].set_ylabel('# S-phase\ncells')

    # sixth row
    sns.barplot(x=df.index, y='num_cells_g', data=df, ax=ax[5,0], color='tab:blue')
    ax[5,0].set_ylabel('# G1/2-phase\ncells')

    # remove the subplot borders from the first 4 rows
    # also remove y-axis ticks and tick labels
    for i in range(4):
        ax[i,0].set_yticklabels([])
        ax[i,0].set_yticks([])
        ax[i,0].spines['top'].set_visible(False)
        ax[i,0].spines['bottom'].set_visible(False)
        ax[i,0].spines['left'].set_visible(False)
        ax[i,0].spines['right'].set_visible(False)

    # remove the xtick labels from all rows
    for i in range(6):
        ax[i,0].set_xticklabels([])
        ax[i,0].set_xticks([])
    
    # plot a titile on the first row
    ax[0,0].set_title('Metacohort Overview')


def plot_clone_embeddings(ax, clone_embeddings):
    pc1_label = 'PC1'
    pc2_label = 'PC2'
    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC2', hue='type', size='num_cells_s', alpha=0.5, ax=ax[6,0], palette=cell_type_cdict)
    ax[6,0].set_xlabel(pc1_label)
    ax[6,0].set_ylabel(pc2_label)
    # only show the 8th and last elements in the legend
    handles, labels = ax[6,0].get_legend_handles_labels()
    handles = [handles[8], handles[-1]]
    labels = [labels[8], labels[-1]]
    ax[6,0].legend(handles, labels, loc='upper right', title='# S cells')
    ax[6,0].set_title('Clone RT PCA embedding\ncolored by type')

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC2', hue='ploidy', size='num_cells_s', alpha=0.5, ax=ax[6,1], palette=ploidy_cdict)
    ax[6,1].set_xlabel(pc1_label)
    ax[6,1].set_ylabel(pc2_label)
    ax[6,1].get_legend().remove()
    ax[6,1].set_title('Clone RT PCA embedding\ncolored by ploidy')

    sns.scatterplot(data=clone_embeddings, x='PC1', y='PC2', hue='signature', size='num_cells_s', alpha=0.5, ax=ax[6,2], palette=signature_cdict)
    ax[6,2].set_xlabel(pc1_label)
    ax[6,2].set_ylabel(pc2_label)
    ax[6,2].get_legend().remove()
    ax[6,2].set_title('Clone RT PCA embedding\ncolored by signature')


def compute_median_and_ci(df, col_prefix_list):
    ''' For each rt profile, compute the median and 95% CI of the posterior distribution. '''
    # initialize empty dataframe
    out_df = pd.DataFrame()
    for col_prefix in col_prefix_list:
        # subset df to just the columns with this prefix
        col_names = [col for col in df.columns if (col.startswith(col_prefix)) or (col in ['chr', 'start', 'end'])]
        temp_df = df[col_names]
        # set chr, start, end as the index
        temp_df = temp_df.set_index(['chr', 'start', 'end'])
        # compute the median and 95% CI of the posterior distribution across all rows in temp_df
        median = temp_df.median(axis=1)
        ci_lower = temp_df.quantile(0.025, axis=1)
        ci_upper = temp_df.quantile(0.975, axis=1)
        median_colname = col_prefix + '_median'
        ci_lower_colname = col_prefix + '_ci_lower'
        ci_upper_colname = col_prefix + '_ci_upper'
        temp_summary_df = pd.DataFrame({median_colname: median, ci_lower_colname: ci_lower, ci_upper_colname: ci_upper}, index=temp_df.index)

        # if out_df is empty, set it equal to temp_summary_df
        if out_df.empty:
            out_df = temp_summary_df
        # otherwise, merge temp_summary_df with out_df
        else:
            out_df = out_df.merge(temp_summary_df, left_index=True, right_index=True)
        
    # reset index so that chr, start, end are now columns
    out_df = out_df.reset_index()

    return out_df


def plot_rt_profile_posteriors(ax, rt_profiles):

    # list of columns to plot
    col_prefix_list = ['global_rt', 'ct_hgsoc', 'ct_tnbc', 'ct_htert', 'sig_fbi', 'wgd_rt']
    subplot_indices = [(7,0), (7,2), (8,0), (8,2), (9,0), (9,2)]
    
    # share the y-axes across all subplots in the subplot indices
    for i,j in subplot_indices:
        ax[i,j].get_shared_y_axes().join(ax[i,j], ax[7,0])

    # plot the RT profiles with the colors corresponding to the column prefix
    for (i,j), col_prefix in zip(subplot_indices, col_prefix_list):
        if col_prefix.startswith('global'):
            color = 'C0'
        elif col_prefix.startswith('ct'):
            color = 'C1'
        elif col_prefix.startswith('sig'):
            color = 'C2'
        elif col_prefix.startswith('wgd') or col_prefix.startswith('ngd'):
            color = 'C3'
        median_colname = col_prefix + '_median'
        ci_lower_colname = col_prefix + '_ci_lower'
        ci_upper_colname = col_prefix + '_ci_upper'
        plot_cell_cn_profile2(ax[i,j], rt_profiles, median_colname, cn_field_name=None, max_cn=None,
                            chromosome=None, s=1, squashy=False, color=color, alpha=1,
                            lines=True, label=None, scale_data=False, rawy=True,
                            rasterized=True, min_ci_field_name=ci_lower_colname, max_ci_field_name=ci_upper_colname)
        title = col_prefix.replace('global_rt', 'Global RT').replace('ct_', 'cell type = ').replace('sig_', 'signature = ').replace('wgd_rt', 'WGD').replace('ngd_rt', 'NGD')
        title = title.replace('htert', 'hTERT').replace('hgsoc', 'HGSOC').replace('tnbc', 'TNBC').replace('fbi', 'FBI')
        ax[i,j].set_title(title)


def plot_main_figure(df, clone_embeddings, feature_coef_df, rt_profiles):
    # create a figure with 10 rows and 4 column of size 12x4
    # have the first 4 rows be half the size as the last 2 rows
    fig, ax = plt.subplots(10, 4, figsize=(12, 10), tight_layout=False, gridspec_kw={'height_ratios': [0.2, 0.2, 0.2, 0.2, 1, 1, 4, 1.5, 1.5, 1.5]})

    # set buffer between rows and columns
    # plt.subplots_adjust(hspace=1.5, wspace=0.5)

    # for the first 6 rows, merge each row into a single axis
    for i in range(6):
        gs = ax[i, 0].get_gridspec()
        for axi in ax[i, :]:
            axi.remove()
        ax[i, 0] = fig.add_subplot(gs[i, :])


    # for the final 3 rows, merge each row into two axes
    for i in range(7, 10):
        gs = ax[i, 0].get_gridspec()
        for axi in ax[i, :]:
            axi.remove()
        ax[i, 0] = fig.add_subplot(gs[i, 0:2])
        ax[i, 2] = fig.add_subplot(gs[i, 2:4])


    # plot the metacohort overview in the first 6 rows
    plot_metacohort_overview(ax, df)

    # plot the clone pca embeddings in the first 3 columns of the 6th row
    plot_clone_embeddings(ax, clone_embeddings)

    # make a barplot of the PCA regression coefficients for each feature
    sns.barplot(data=feature_coef_df, x='metafeature', y='log10_abs_coef', ax=ax[6,3], palette=feature_cdict)
    ax[6,3].set_ylabel('log10(abs(Coefficient))')
    ax[6,3].set_xlabel('Feature')
    ax[6,3].set_title('Multivariate linear regression of PCs\nusing K-fold CV')

    # plot the rt profile posteriors in the final 3 columns of the 7th-9th rows
    plot_rt_profile_posteriors(ax, rt_profiles)

    return fig


def main():
    argv = get_args()

    # load in the data
    df = pd.read_csv(argv.features)
    table = pd.read_csv(argv.metacohort, sep='\t')
    clone_embeddings = pd.read_csv(argv.pca_embeddings)
    feature_coef_df = pd.read_csv(argv.pca_regression)
    rt_profiles = pd.read_csv(argv.rt_profiles)

    # drop all the one-hot-encoded columns from df that begin with type_ or signature_ 
    df = df.loc[:,~df.columns.str.startswith('type_')]
    df = df.loc[:,~df.columns.str.startswith('signature_')]

    # drop the cn_path and rt_path columns from table
    table = table.drop(['cn_path', 'rt_path', 'condition'], axis=1)

    # merge the two dataframes on the dataset column
    df = df.merge(table, on='dataset')

    # reset ploidy to be an integer
    clone_embeddings['ploidy'] = clone_embeddings['ploidy'].astype(int)

    # list of all the rt profiles that need to be processed
    col_prefix_list = ['global_rt', 'ct_hgsoc', 'ct_tnbc', 'ct_htert', 'ct_ov2295', 'ct_t47d', 'ct_gm18507', 'sig_fbi', 'sig_hrd', 'sig_td', 'wgd_rt', 'ngd_rt']

    # compute the median and 95% CI of the posterior distribution for each rt profile
    rt_profiles = compute_median_and_ci(rt_profiles, col_prefix_list)

    # set the chromosome column as a categorical variable
    rt_profiles['chr'] = rt_profiles['chr'].astype('str').astype('category')

    # plot the main figure
    fig = plot_main_figure(df, clone_embeddings, feature_coef_df, rt_profiles)

    # save the figure
    fig.savefig(argv.output, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()

