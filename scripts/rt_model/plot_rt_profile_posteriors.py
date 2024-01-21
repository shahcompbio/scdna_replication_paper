import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser
from scdna_replication_tools.plot_utils import plot_cell_cn_profile2


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='table of rt profile posteriors')
    p.add_argument('--noX', action='store_true', help='binary indicator to include for version of the model where chrX is excluded')
    p.add_argument('output_csv', type=str, help='csv of mean, median, upper & lowe CI intervals from posterior rt profiles')
    p.add_argument('output_pdf', type=str, help='plot of rt profiles sampled from the posterior distribution')

    return p.parse_args()


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


def main():
    argv = get_args()

    # read in the posterior rt profiles
    rt_profiles = pd.read_csv(argv.input)

    # list of all the rt profiles that need to be processed
    col_prefix_list = ['global_rt', 'ct_hgsoc', 'ct_tnbc', 'ct_htert', 'ct_ov2295', 'ct_t47d', 'ct_gm18507', 'sig_fbi', 'sig_hrd', 'sig_td', 'wgd_rt', 'ngd_rt']

    # compute the median and 95% CI of the posterior distribution for each rt profile
    rt_profiles = compute_median_and_ci(rt_profiles, col_prefix_list)

    # set the chromosome column as a categorical variable
    rt_profiles['chr'] = rt_profiles['chr'].astype('str').astype('category')

    # plot the genome-wide learned parameters
    fig, ax = plt.subplots(6, 2, figsize=(12, 8), tight_layout=True, sharey=True)
    ax = ax.flatten()

    for i, col_prefix in enumerate(col_prefix_list):
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
        plot_cell_cn_profile2(ax[i], rt_profiles, median_colname, cn_field_name=None, max_cn=None,
                            chromosome=None, s=1, squashy=False, color=color, alpha=1,
                            lines=True, label=None, scale_data=False, rawy=True,
                            rasterized=True, min_ci_field_name=ci_lower_colname, max_ci_field_name=ci_upper_colname)
        title = col_prefix.replace('global_rt', 'Global RT').replace('ct_', 'cell type = ').replace('sig_', 'signature = ').replace('wgd_rt', 'WGD').replace('ngd_rt', 'NGD')
        if argv.noX:
            title = title + ' (no chrX)'
        ax[i].set_title(title)

    fig.savefig(argv.output_pdf, dpi=300, bbox_inches='tight')
    rt_profiles.to_csv(argv.output_csv, index=False)


if __name__ == '__main__':
    main()
