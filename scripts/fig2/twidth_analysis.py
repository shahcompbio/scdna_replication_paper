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

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells with true and inferred scRT data')
    p.add_argument('dataset')
    p.add_argument('A', type=float, help='steepness of inflection point when drawing RT state')
    p.add_argument('nb_r', type=float, help='amount of sequencing noise')
    p.add_argument('rt_col', help='true replication time of each bin (i.e. mcf7rt')
    p.add_argument('frac_rt_col', help='inferred fraction replicated for each cell')
    p.add_argument('true_frac_col', help='true fraction replicated for each cell')
    p.add_argument('rep_state', help='inferred replication state for each bin')
    p.add_argument('true_rep_state', help='true replication state for each bin')
    p.add_argument('infer_mode', help='pyro model or bulk')
    p.add_argument('output_heatmap', help='heatmap comparing true and inferred rt_state values with T-width superimposed')
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


# helper functions for plotting heatmaps
def plot_colorbar(ax, color_mat, title=None):
    ax.imshow(np.array(color_mat)[::-1, np.newaxis], aspect='auto', origin='lower')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    if title is not None:
        ax.set_title(title)


def plot_color_legend(ax, color_map, title=None):
    legend_elements = []
    for name, color in color_map.items():
        legend_elements.append(Patch(facecolor=color, label=name))
    ax.legend(handles=legend_elements, loc='center left', title=title)
    ax.grid(False)
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])


def make_color_mat_float(values, palette_color):
    """
    Make a color_mat for a 0-1 float array `values` and a
    corresponding color pallete.
    """
    pal = plt.get_cmap(palette_color)
    color_mat = []
    for val in values:
        color_mat.append(pal(val))
    color_dict = {0: pal(0.0), 1: pal(1.0)}
    return color_mat, color_dict


def get_rt_cmap():
    rt_colors = {0: '#552583', 1: '#FDB927'}
    color_list = []
    for i in [0, 1]:
        color_list.append(rt_colors[i])
    return ListedColormap(color_list)


def plot_true_vs_inferred_rt_state(df, true_Tw, inferred_Tw, title_second_line, argv):
    df = df.copy()

    # create mapping of clones
    cluster_col = 'cluster_id'
    clone_col = 'clone_id'
    clone_dict = dict([(y,x+1) for x,y in enumerate(sorted(df[clone_col].unique()))])
    df[cluster_col] = df[clone_col]
    df = df.replace({cluster_col: clone_dict})

    secondary_sort_column = argv.true_frac_col
    secondary_sort_label = 'Frac Rep'

    fig = plt.figure(figsize=(14, 7))
    
    ax0 = fig.add_axes([0.12,0.0,0.38,1.])
    rt_cmap = get_rt_cmap()
    plot_data0 = plot_clustered_cell_cn_matrix(ax0, df, argv.true_rep_state, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column, cmap=rt_cmap)
    ax0.set_title('True scRT, T-width: {}\n{}'.format(round(true_Tw, 3), title_second_line))

    ax1 = fig.add_axes([0.62,0.0,0.38,1.])
    plot_data1 = plot_clustered_cell_cn_matrix(ax1, df, argv.rep_state, cluster_field_name=cluster_col, secondary_field_name=secondary_sort_column, cmap=rt_cmap)
    ax1.set_title('Inferred scRT, T-width: {}\n{}'.format(round(inferred_Tw, 3), title_second_line))

    if len(clone_dict) > 1:
        # annotate the clones for G1-phase cells
        cell_ids = plot_data0.columns.get_level_values(0).values
        cluster_ids0 = plot_data0.columns.get_level_values(1).values
        color_mat0, color_map0 = cncluster.get_cluster_colors(cluster_ids0, return_map=True)

        # get list of color pigments in the same order as clone_dict
        colors_used0 = []
        for c in color_mat0:
            if c not in colors_used0:
                colors_used0.append(c)

        # match clone IDs to color pigments
        clones_to_colors0 = {}
        for i, key in enumerate(clone_dict.keys()):
            clones_to_colors0[key] = colors_used0[i]

        # get array of secondary_sort_column values that that match the cell_id order
        condensed_cn = df[['cell_id', secondary_sort_column]].drop_duplicates()
        secondary_array = []
        for cell in cell_ids:
            s = condensed_cn[condensed_cn['cell_id'] == cell][secondary_sort_column].values[0]
            secondary_array.append(s)

        # make color mat according to secondary array
        secondary_color_mat, secondary_to_colors = make_color_mat_float(secondary_array, 'Blues')

        # create color bar that shows clone id for each row in heatmap
        ax = fig.add_axes([0.09,0.0,0.03,1.])
        plot_colorbar(ax, color_mat0)

        # create color bar that shows secondary sort value for each row in heatmap
        ax = fig.add_axes([0.06,0.0,0.03,1.])
        plot_colorbar(ax, secondary_color_mat)

        # create legend to match colors to clone ids
        ax = fig.add_axes([0.0,0.75,0.04,0.25])
        plot_color_legend(ax, clones_to_colors0, title='Clone ID')

        # create legend to match colors to secondary sort values
        ax = fig.add_axes([0.0,0.5,0.04,0.25])
        plot_color_legend(ax, secondary_to_colors, title=secondary_sort_label)

        # annotate the clones for S-phase cells.. using the same colors as G1 clones
        cluster_ids1 = plot_data1.columns.get_level_values(1).values
        color_mat1 = cncluster.get_cluster_colors(cluster_ids1, color_map=color_map0)

        # match clone IDs to color pigments
        clones_to_colors1 = {}
        for i, key in enumerate(clone_dict.keys()):
            clones_to_colors1[key] = colors_used0[i]

        # create color bar that shows clone id for each row in heatmap
        ax = fig.add_axes([0.59,0.0,0.03,1.])
        plot_colorbar(ax, color_mat1)

        # create color bar that shows secondary sort value for each row in heatmap
        ax = fig.add_axes([0.56,0.0,0.03,1.])
        plot_colorbar(ax, secondary_color_mat)

        # create legend to match colors to clone ids
        ax = fig.add_axes([0.5,0.75,0.04,0.25])
        plot_color_legend(ax, clones_to_colors1, title='Clone ID')

         # create legend to match colors to secondary sort values
        ax = fig.add_axes([0.5,0.5,0.04,0.25])
        plot_color_legend(ax, secondary_to_colors, title=secondary_sort_label)

    return fig



def main():
    argv = get_args()
    df = pd.read_csv(argv.cn_s, sep='\t')

    # set chr column to category
    df.chr = df.chr.astype('str')
    df.chr = df.chr.astype('category')

    # compute time from scheduled replication for each bin
    df['rt_hours'] = ((df[argv.rt_col] / 10.0) - 10.) * -1.
    df['time_from_scheduled_rt'] = df['rt_hours'] - (df[argv.frac_rt_col] * 10.0)
    df['true_time_from_scheduled_rt'] = df['rt_hours'] - (df[argv.true_frac_col] * 10.0)

    fig, ax = plt.subplots(1, 2, figsize=(10, 5), tight_layout=True)
    ax = ax.flatten()

    title_second_line = 'dataset: {}, A: {}, nb_r: {}, infer: {}'.format(
        argv.dataset, argv.A, argv.nb_r, argv.infer_mode
    )
    Tw = compute_and_plot_twidth(df, column='time_from_scheduled_rt', rep_state=argv.rep_state, title='Inferred scRT heterogeneity\n{}'.format(title_second_line), ax=ax[1])
    true_Tw = compute_and_plot_twidth(df, column='true_time_from_scheduled_rt', rep_state=argv.true_rep_state, title='True scRT heterogeneity\n{}'.format(title_second_line), ax=ax[0])

    fig.savefig(argv.output_curves, bbox_inches='tight')

    fig = plot_true_vs_inferred_rt_state(df, true_Tw, Tw, title_second_line, argv)
    
    fig.savefig(argv.output_heatmap, bbox_inches='tight')


if __name__ == '__main__':
    main()
