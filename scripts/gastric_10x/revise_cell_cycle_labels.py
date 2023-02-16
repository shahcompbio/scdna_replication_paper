import pandas as pd
import numpy as np
import statsmodels.api as sm
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('input_s', type=str, help='pyro model output for s-phase cells')
    p.add_argument('input_g', type=str, help='pyro model output for g1/2-phase cells')
    p.add_argument('frac_rt_col', type=str, help='column name for the fraction of replicated bins in a cell')
    p.add_argument('rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('cn_col', type=str, help='column name for inferred cn state of each bin in pyro model')
    p.add_argument('rpm_col', type=str, help='column name for reads per million')
    p.add_argument('output_s', type=str, help='pyro model for s-phase cells after filtering')
    p.add_argument('output_g', type=str, help='pyro model results for nonreplicating cells that get removed')
    p.add_argument('output_lowqual', type=str, help='pyro model results for low quality cells that get removed')

    return p.parse_args()


def autocorr(data, min_lag=10, max_lag=50):
    acorr = sm.tsa.acf(data, nlags=max_lag)
    return np.mean(acorr[min_lag-1:])


def breakpoints(data):
    temp_diff = np.diff(data)
    return sum(np.where(temp_diff!=0, 1, 0))


def compute_cell_frac(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state'):
    ''' Compute the fraction of replicated bins for all cells in `cn`, noting which cells have extreme values. '''
    cn['extreme_cell_frac'] = False
    for cell_id, cell_cn in cn.groupby('cell_id'):
        temp_rep = cell_cn[rep_state_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[cell_cn.index, frac_rt_col] = temp_frac

        if temp_frac > 0.95 or temp_frac < 0.05:
            cn.loc[cell_cn.index, 'extreme_cell_frac'] = True
    return cn


def remove_nonreplicating_cells(cn, frac_rt_col='cell_frac_rep'):
    # use extreme_cell_frac status and cell cycle classifier features to nominate "bad" non-replicating cells
    extreme_cells = cn.loc[(cn['extreme_cell_frac']==True)]
    # bad_cells_df = extreme_cells.loc[(extreme_cells['corrected_breakpoints']<0.0) | (extreme_cells['corrected_madn']<0.0)]
    # bad_cells = bad_cells_df.cell_id.unique()
    bad_cells = extreme_cells.cell_id.unique()

    cn_good = cn[~cn['cell_id'].isin(bad_cells)].reset_index(drop=True)
    cn_bad = cn[cn['cell_id'].isin(bad_cells)].reset_index(drop=True)

    return cn_good, cn_bad


def compute_quality_features(cn, rep_state_col='model_rep_state', cn_state_col='model_cn_state', rpm_col='rpm'):
    cell_metrics = []
    for cell_id, cell_cn in cn.groupby('cell_id'):
        # compute read depth autocorrelation
        rpm_auto = autocorr(cell_cn[rpm_col].values)
        # compute replication state autocorrelation
        rep_auto = autocorr(cell_cn[rep_state_col].values)
        # compute number of breakpoints for inferred CN
        cn_bk = breakpoints(cell_cn[cn_state_col].values)
        # compute number of breakpoints for inferred rep states
        rep_bk = breakpoints(cell_cn[rep_state_col].values)
        # compute fraction of genome with inferred CN=0
        frac_cn0 = sum(cell_cn[cn_state_col]==0) / cell_cn.shape[0]

        temp_df = pd.DataFrame({
            'cell_id': [cell_id], 'rpm_auto': [rpm_auto], 'rep_auto': [rep_auto],
            'cn_bk': [cn_bk], 'rep_bk': [rep_bk], 'frac_cn0': [frac_cn0]
        })
        cell_metrics.append(temp_df)

    cell_metrics = pd.concat(cell_metrics, ignore_index=True)
    
    # normalize autocorrelations by mean values
    rpm_auto_mean = np.mean(cell_metrics['rpm_auto'].values)
    cell_metrics['rpm_auto_norm'] = cell_metrics['rpm_auto'] - rpm_auto_mean
    rep_auto_mean = np.mean(cell_metrics['rep_auto'].values)
    cell_metrics['rep_auto_norm'] = cell_metrics['rep_auto'] - rep_auto_mean

    # merge quality features into cn
    cn = pd.merge(cn, cell_metrics)
    
    return cn


def remove_low_quality_cells(cn):
    # set thresholds on which cells are high or low quality
    # note that this thresholding takes place after removing cells with extreme fractions of replicated bins
    low_qual_cells = cn.loc[(cn['rpm_auto']>0.5) | (cn['rep_auto']>0.2) | (cn['frac_cn0']>0.05)].cell_id.unique()

    cn_good = cn[~cn['cell_id'].isin(low_qual_cells)].reset_index(drop=True)
    cn_bad = cn[cn['cell_id'].isin(low_qual_cells)].reset_index(drop=True)

    return cn_good, cn_bad


if __name__ == '__main__':
    argv = get_args()

    # load in unfiltered S-phase cells
    cn_s = pd.read_csv(argv.input_s, sep='\t')
    cn_g = pd.read_csv(argv.input_g, sep='\t')
    cn = pd.concat([cn_s, cn_g], ignore_index=True)

    # compute the fraction of replicated bins within each cell
    cn = compute_cell_frac(cn, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.rep_col)

    # compute autocorrelation features to see which cells are truly low quality
    cn = compute_quality_features(cn, rep_state_col=argv.rep_col, cn_state_col=argv.cn_col, rpm_col=argv.rpm_col)

    # remove set of cells that appear to be nonreplicating
    cn_filtered, cn_nonrep = remove_nonreplicating_cells(cn, frac_rt_col=argv.frac_rt_col)

    # remove set of cells that appear to be low quality
    cn_filtered, cn_lowqual = remove_low_quality_cells(cn_filtered)

    # return the filtered set of S-phase cells
    cn_filtered.to_csv(argv.output_s, sep='\t', index=False)

    # return the set of cells that don't appear to be replicating
    cn_nonrep.to_csv(argv.output_g, sep='\t', index=False)

    # return the cells that look to be low quality
    cn_lowqual.to_csv(argv.output_lowqual, sep='\t', index=False)

