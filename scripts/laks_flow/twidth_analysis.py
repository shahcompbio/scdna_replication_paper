from argparse import ArgumentParser
import pandas as pd
import matplotlib.pyplot as plt
from scdna_replication_tools.calculate_twidth import compute_and_plot_twidth


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_T47D', help='input long-form copy number dataframe for T47D S-phase cells with inferred scRT')
    p.add_argument('cn_GM18507', help='input long-form copy number dataframe for GM18507 S-phase cells with inferred scRT')
    p.add_argument('cn_all', help='input long-form copy number dataframe for T47D + GM18507 S-phase cells with inferred scRT')
    p.add_argument('pseduobulk', help='RT pseudobulk for this dataset')
    p.add_argument('frac_rt_col', help='inferred fraction replicated for each cell')
    p.add_argument('rep_state', help='inferred replication state for each bin')
    p.add_argument('output_tsv', help='table of all the computed t_width values')
    p.add_argument('output_curves', help='T-width curves of true and inferred rt_states')

    return p.parse_args()


def compute_loci_frac(cn, rt_state='model_rep_state', col_name='loci_frac_rep'):
    ''' Compute the fraction of replicated bins at each locus '''
    for (chrom, start), loci_cn in cn.groupby(['chr', 'start']):
        temp_rep = loci_cn[rt_state].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[loci_cn.index, col_name] = temp_frac
    return cn


def run_twidth_analysis(df, rt_col, frac_rt_col, rep_state, title_second_line='', per_cell=False, alpha=1, ax=None):
    df = df.copy()
    
    # set chr column to category
    df['chr'] = df['chr'].astype('str')
    df['chr'] = df['chr'].astype('category')

    # compute time from scheduled replication for each bin
    df['rt_hours'] = ((df[rt_col] * 10.0) - 10.) * -1.
    df['time_from_scheduled_rt'] = df['rt_hours'] - (df[frac_rt_col] * 10.0)

    Tw = compute_and_plot_twidth(
        df, tfs_col='time_from_scheduled_rt', rs_col=rep_state,
        title='Cellular scRT heterogeneity\n{}'.format(title_second_line),
        per_cell=per_cell, alpha=alpha, ax=ax
    )
    
    return Tw


def main():
    argv = get_args()
    
    # load the data
    cn_t = pd.read_csv(argv.cn_T47D, sep='\t')
    cn_gm = pd.read_csv(argv.cn_GM18507, sep='\t')
    cn_all = pd.read_csv(argv.cn_all, sep='\t')
    rt_bulks = pd.read_csv(argv.pseduobulk, sep='\t')

    # merge rt bulk info into the main cn dataframes
    cn_all = pd.merge(cn_all, rt_bulks)
    cn_t = pd.merge(cn_t, rt_bulks)
    cn_gm = pd.merge(cn_gm, rt_bulks)
    cn_split = pd.concat([cn_t, cn_gm], ignore_index=True)

    # compute a pseudobulk RT column that is the average of the two cell lines
    cn_all = compute_loci_frac(cn_all, rt_state=argv.rep_state, col_name='rt_merged_all')
    cn_split = compute_loci_frac(cn_split, rt_state=argv.rep_state, col_name='rt_split_all')

    t_width_df = []

    # generate one subplot per cell line, model condition and per_cell status
    fig, ax = plt.subplots(4, 3, figsize=(12,16), tight_layout=True)
    ax = ax.flatten()

    # plots for per_cell==False
    ax[0], Tw = run_twidth_analysis(cn_gm, 'rt_split_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 Split', ax=ax[0])
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['split'], 'per_cell': [False], 'T-width': [Tw]}))

    ax[1], Tw = run_twidth_analysis(cn_t, 'rt_split_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D Split', ax=ax[1])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['split'], 'per_cell': [False], 'T-width': [Tw]}))

    ax[2], Tw = run_twidth_analysis(cn_split, 'rt_split_all', argv.frac_rt_col, argv.rep_state, title_second_line='Split T47D + GM18507', ax=ax[2])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['split'], 'per_cell': [False], 'T-width': [Tw]}))

    ax[3], Tw = run_twidth_analysis(cn_gm, 'rt_merged_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 merged', ax=ax[3])
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['merged'], 'per_cell': [False], 'T-width': [Tw]}))

    ax[4], Tw = run_twidth_analysis(cn_t, 'rt_merged_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D merged', ax=ax[4])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['merged'], 'per_cell': [False], 'T-width': [Tw]}))

    ax[5], Tw = run_twidth_analysis(cn_all, 'rt_merged_all', argv.frac_rt_col, argv.rep_state, title_second_line='merged T47D + GM18507', ax=ax[5])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['merged'], 'per_cell': [False], 'T-width': [Tw]}))

    # plots for per-cell==True
    ax[6], Tw = run_twidth_analysis(cn_gm, 'rt_split_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 Split', ax=ax[6], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['split'], 'per_cell': [True], 'T-width': [Tw]}))

    ax[7], Tw = run_twidth_analysis(cn_t, 'rt_split_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D Split', ax=ax[7], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['split'], 'per_cell': [True], 'T-width': [Tw]}))

    ax[8], Tw = run_twidth_analysis(cn_split, 'rt_split_all', argv.frac_rt_col, argv.rep_state, title_second_line='Split T47D + GM18507', ax=ax[8], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['split'], 'per_cell': [True], 'T-width': [Tw]}))

    ax[9], Tw = run_twidth_analysis(cn_gm, 'rt_merged_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 merged', ax=ax[9], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['merged'], 'per_cell': [True], 'T-width': [Tw]}))

    ax[10], Tw = run_twidth_analysis(cn_t, 'rt_merged_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D merged', ax=ax[10], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['merged'], 'per_cell': [True], 'T-width': [Tw]}))

    ax[11], Tw = run_twidth_analysis(cn_all, 'rt_merged_all', argv.frac_rt_col, argv.rep_state, title_second_line='merged T47D + GM18507', ax=ax[11], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['merged'], 'per_cell': [True], 'T-width': [Tw]}))

    # save figure of all the t-width curves
    fig.savefig(argv.output_curves, bbox_inches='tight')

    # save a table of all the computed T-width values
    t_width_df = pd.concat(t_width_df, ignore_index=True)
    t_width_df.to_csv(argv.output_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
