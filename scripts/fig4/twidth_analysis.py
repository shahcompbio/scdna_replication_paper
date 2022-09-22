from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from matplotlib.colors import ListedColormap
from scipy.optimize import curve_fit
from scgenome import cncluster
from matplotlib.patches import Patch


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


def calc_pct_replicated_per_time_bin(cn, column='time_from_scheduled_rt', per_cell=False, query2=None, rep_state='rt_state'):
    intervals = np.linspace(-10, 10, 201)
    time_bins = []
    pct_reps = []
    for i in range(200):
        a = intervals[i]
        b = intervals[i+1]
        temp_cn = cn.query("{col} < {b} & {col} >= {a}".format(a=a, b=b, col=column))
        if query2:
            temp_cn = temp_cn.query(query2)
        if temp_cn.shape[0] > 0:
            if per_cell:
                for cell_id, chunk_cn in temp_cn.groupby('cell_id'):
                    if chunk_cn.shape[0] > 0:
                        percent_replicated = sum(chunk_cn[rep_state]) / len(chunk_cn[rep_state])
                        pct_reps.append(percent_replicated)
                        time_bins.append(a)
            else:
                percent_replicated = sum(temp_cn[rep_state]) / len(temp_cn[rep_state])
                pct_reps.append(percent_replicated)
                time_bins.append(a)
    return time_bins, pct_reps


# helper functions for computing T width
def sigmoid(x, x0, k, b):
    y = 1 / (1 + np.exp(-k*(x-x0)))+b
    return y

def inv_sigmoid(y, x0, k, b):
    temp = (1 / (y-b)) - 1
    x = (np.log(temp) / -k) + x0
    return x

def fit_sigmoid(xdata, ydata):
    p0 = [np.median(xdata), 1, 0.0] # this is an mandatory initial guess
    popt, pcov = curve_fit(sigmoid, xdata, ydata, p0, method='dogbox')
    return popt, pcov

def calc_t_width(popt, low=0.25, high=0.75):
    right_time = inv_sigmoid(low, *popt)
    left_time = inv_sigmoid(high, *popt)
    t_width = right_time - left_time
    return t_width, left_time, right_time

# helper functions to use if you want to calculate T-width via linear instead of sigmoid regression
def linear(x, m, b):
    y = m * x + b
    return y

def inv_linear(y, m, b):
    x = (y - b) / m
    return x

def fit_linear(xdata, ydata):
    p0 = [-1.0, -1.0] # this is an mandatory initial guess
    popt, pcov = curve_fit(linear, xdata, ydata, p0)
    return popt, pcov

def calc_linear_t_width(popt, low=0.25, high=0.75):
    right_time = inv_linear(low, *popt)
    left_time = inv_linear(high, *popt)
    t_width = right_time - left_time
    return t_width, left_time, right_time


# superimpose T-width lines onto sigmoid function & data
def plot_cell_variability(xdata, ydata, popt=None, left_time=None, right_time=None, t_width=None,
                          alpha=1, title='Cell-to-cell variabilty', curve='sigmoid', ax=None):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    
    ax.scatter(xdata, ydata, label='data', alpha=alpha)
    if popt is not None:
        x = np.linspace(-10, 10, 1000)
        if curve == 'sigmoid':
            y = sigmoid(x, *popt)
        elif curve == 'linear':
            y = linear(x, *popt)
        ax.plot(x, y, color='r', label='fit')
        ax.axhline(y=0.75, color='k', linestyle='--')
        ax.axhline(y=0.25, color='k', linestyle='--')
        ax.axvline(x=left_time, color='k', linestyle='--')
        ax.axvline(x=right_time, color='k', linestyle='--', label='T_width={}'.format(round(t_width, 3)))
    ax.set_xlabel('time from scheduled replication (h)')
    ax.set_ylabel('% replicated')
    ax.set_title(title)
    ax.legend(loc='best')


def compute_and_plot_twidth(cn, column='time_from_scheduled_rt', per_cell=False, query2=None, rep_state='rt_state',
                            alpha=1, title='Cell-to-cell variabilty', curve='sigmoid', ax=None):
    time_bins, pct_reps = calc_pct_replicated_per_time_bin(cn, per_cell=per_cell, column=column, query2=query2, rep_state=rep_state)
    if curve == 'sigmoid':
        popt, pcov = fit_sigmoid(time_bins, pct_reps)
        t_width, right_time, left_time = calc_t_width(popt)
    elif curve == 'linear':
        popt, pcov = fit_linear(time_bins, pct_reps)
        t_width, right_time, left_time = calc_linear_t_width(popt)
    plot_cell_variability(time_bins, pct_reps, popt,
                          left_time, right_time, t_width,
                          alpha=alpha, title=title, curve=curve, ax=ax)
    
    return t_width


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
        df, column='time_from_scheduled_rt', rep_state=rep_state,
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
    cn_all = compute_loci_frac(cn_all, rt_state=argv.rep_state, col_name='rt_joint_all')
    cn_split = compute_loci_frac(cn_split, rt_state=argv.rep_state, col_name='rt_split_all')

    t_width_df = []

    # generate one subplot per cell line, model condition and per_cell status
    fig, ax = plt.subplots(2, 3, figsize=(12,8), tight_layout=True)
    ax = ax.flatten()

    # plots for per_cell==False
    Tw = run_twidth_analysis(cn_gm, 'rt_split_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 Split', ax=ax[0])
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['split'], 'per_cell': [False], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_t, 'rt_split_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D Split', ax=ax[1])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['split'], 'per_cell': [False], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_split, 'rt_split_all', argv.frac_rt_col, argv.rep_state, title_second_line='Split T47D + GM18507', ax=ax[2])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['split'], 'per_cell': [False], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_gm, 'rt_joint_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 Joint', ax=ax[3])
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['joint'], 'per_cell': [False], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_t, 'rt_joint_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D Joint', ax=ax[4])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['joint'], 'per_cell': [False], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_split, 'rt_joint_all', argv.frac_rt_col, argv.rep_state, title_second_line='Joint T47D + GM18507', ax=ax[5])
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['joint'], 'per_cell': [False], 'T-width': [Tw]}))

    # plots for per-cell==True
    Tw = run_twidth_analysis(cn_gm, 'rt_split_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 Split', ax=ax[6], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['split'], 'per_cell': [True], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_t, 'rt_split_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D Split', ax=ax[7], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['split'], 'per_cell': [True], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_split, 'rt_split_all', argv.frac_rt_col, argv.rep_state, title_second_line='Split T47D + GM18507', ax=ax[8], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['split'], 'per_cell': [True], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_gm, 'rt_joint_GM18507', argv.frac_rt_col, argv.rep_state, title_second_line='GM18507 Joint', ax=ax[9], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['GM18507'], 'model_condition': ['joint'], 'per_cell': [True], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_t, 'rt_joint_T47D', argv.frac_rt_col, argv.rep_state, title_second_line='T47D Joint', ax=ax[10], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D'], 'model_condition': ['joint'], 'per_cell': [True], 'T-width': [Tw]}))

    Tw = run_twidth_analysis(cn_split, 'rt_joint_all', argv.frac_rt_col, argv.rep_state, title_second_line='Joint T47D + GM18507', ax=ax[11], per_cell=True, alpha=0.1)
    t_width_df.append(pd.DataFrame({'cell_line': ['T47D+GM18507'], 'model_condition': ['joint'], 'per_cell': [True], 'T-width': [Tw]}))

    # save figure of all the t-width curves
    fig.savefig(argv.output_curves, bbox_inches='tight')

    # save a table of all the computed T-width values
    t_width_df = pd.concat(t_width_df, ignore_index=True)
    t_width_df.to_csv(argv.output_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
