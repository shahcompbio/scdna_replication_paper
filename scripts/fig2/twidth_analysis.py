from argparse import ArgumentParser
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scgenome.cnplot import plot_clustered_cell_cn_matrix
from matplotlib.colors import ListedColormap
from scipy.optimize import curve_fit


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('dataset')
    p.add_argument('sigma1', type=float, help='noise of read depth profiles')
    p.add_argument('gc_slope', type=float, help='slope of linear GC bias')
    p.add_argument('gc_int', type=float, help='intercept of linear GC bias')
    p.add_argument('A', type=float, help='steepness of inflection point when drawing RT state')
    p.add_argument('s_time_dist', help='distribution of S-phase cells captured (normal or uniform)')
    p.add_argument('output_heatmap', help='heatmap comparing true and inferred rt_state values with T-width superimposed')
    p.add_argument('output_curves', help='T-width curves of true and inferred rt_states')

    return p.parse_args()


def get_rt_cmap():
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list)


def calc_pct_replicated_per_time_bin(cn, column='time_from_scheduled_rt', per_cell=False, query2=None):
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
                        percent_replicated = sum(chunk_cn.rt_state) / len(chunk_cn.rt_state)
                        pct_reps.append(percent_replicated)
                        time_bins.append(a)
            else:
                percent_replicated = sum(temp_cn.rt_state) / len(temp_cn.rt_state)
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


def compute_and_plot_twidth(cn, column='time_from_scheduled_rt', per_cell=False, query2=None,
                            alpha=1, title='Cell-to-cell variabilty', curve='sigmoid', ax=None):
    time_bins, pct_reps = calc_pct_replicated_per_time_bin(cn, per_cell=per_cell, column=column, query2=query2)
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


def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # compute time from scheduled replication for each bin
    df['mcf7rt_hours'] = ((df['mcf7rt'] / 10.0) - 10.) * -1.
    df['time_from_scheduled_rt'] = df['mcf7rt_hours'] - (df['frac_rt'] * 10.0)
    df['true_time_from_scheduled_rt'] = df['mcf7rt_hours'] - (df['true_frac_rt'] * 10.0)


    fig, ax = plt.subplots(1, 2, figsize=(10, 5), tight_layout=True)
    ax = ax.flatten()

    title_second_line = 'dataset: {}, A: {}, sigma1: {}\ngc_int: {}, gc_slope: {}, s_time_dist: {}'.format(
        argv.dataset, argv.A, argv.sigma1, argv.gc_int, argv.gc_slope, argv.s_time_dist
    )
    Tw = compute_and_plot_twidth(df, column='time_from_scheduled_rt', title='Inferred scRT heterogeneity\n{}'.format(title_second_line), ax=ax[1])
    true_Tw = compute_and_plot_twidth(df, column='true_time_from_scheduled_rt', title='True scRT heterogeneity\n{}'.format(title_second_line), ax=ax[0])

    fig.savefig(argv.output_curves)

    fig, ax = plt.subplots(1, 2, figsize=(14, 7), tight_layout=True)
    ax = ax.flatten()

    rt_cmap = get_rt_cmap()
    plot_clustered_cell_cn_matrix(ax[0], df, 'true_rt_state', cluster_field_name='clone_id', secondary_field_name='true_frac_rt', cmap=rt_cmap)
    plot_clustered_cell_cn_matrix(ax[1], df, 'rt_state', cluster_field_name='clone_id', secondary_field_name='true_frac_rt', cmap=rt_cmap)

    ax[0].set_title('True scRT, T-width: {}\n{}'.format(round(true_Tw, 3), title_second_line))
    ax[1].set_title('Inferred scRT, T-width: {}\n{}'.format(round(Tw, 3), title_second_line))

    fig.savefig(argv.output_heatmap)


if __name__ == '__main__':
    main()
