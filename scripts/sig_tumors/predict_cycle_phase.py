import pandas as pd
# from scdna_replication_tools.predict_cycle_phase import predict_cycle_phase
from scdna_replication_tools.predict_cycle_phase import compute_cell_frac, remove_nonreplicating_cells, compute_quality_features, remove_low_quality_cells
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('input_s', type=str, help='pyro model output for s-phase cells')
    p.add_argument('input_g', type=str, help='pyro model output for g1/2-phase cells')
    p.add_argument('frac_rt_col', type=str, help='column name for the fraction of replicated bins in a cell')
    p.add_argument('rep_state_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('cn_state_col', type=str, help='column name for inferred cn state of each bin in pyro model')
    p.add_argument('rpm_col', type=str, help='column name for reads per million')
    p.add_argument('output_s', type=str, help='pyro model for s-phase cells after filtering')
    p.add_argument('output_g', type=str, help='pyro model results for nonreplicating cells that get removed')
    p.add_argument('output_lowqual', type=str, help='pyro model results for low quality cells that get removed')

    return p.parse_args()



def predict_cycle_phase(cn, frac_rt_col='cell_frac_rep', rep_state_col='model_rep_state', cn_state_col='model_cn_state', rpm_col='rpm'):
    # compute fraction of replicated bins for all cells in `cn`
    cn = compute_cell_frac(cn, frac_rt_col=frac_rt_col, rep_state_col=rep_state_col)

    # compute autocorrelation features to see which cells are truly low quality
    cn = compute_quality_features(cn, rep_state_col=rep_state_col, cn_state_col=cn_state_col, rpm_col=rpm_col)

    # remove set of cells that appear to be nonreplicating
    cn_s_lq, cn_g = remove_nonreplicating_cells(cn, frac_rt_col=frac_rt_col)

    # remove set of cells that appear to be low quality
    cn_s, cn_lq = remove_low_quality_cells(cn_s_lq)

    # add columns to denote the PERT predicted phases
    cn_s['PERT_phase'] = 'S'
    cn_g['PERT_phase'] = 'G1/2'
    cn_lq['PERT_phase'] = 'LQ'

    return cn_s, cn_g, cn_lq



if __name__ == '__main__':
    argv = get_args()

    # load in unfiltered S-phase cells
    cn_s = pd.read_csv(argv.input_s)
    cn_g = pd.read_csv(argv.input_g)

    # concatenate the two dataframes
    cn = pd.concat([cn_s, cn_g], ignore_index=True)

    # convert the 'chr' column to a string and then categorical
    cn.chr = cn.chr.astype(str).astype('category')

    # predict the cycle phase for each cell
    cn_s, cn_g, cn_lq = predict_cycle_phase(
        cn, frac_rt_col=argv.frac_rt_col, rep_state_col=argv.rep_state_col, 
        cn_state_col=argv.cn_state_col, rpm_col=argv.rpm_col
    )

    # return the filtered set of S-phase cells
    cn_s.to_csv(argv.output_s, index=False)

    # return the set of cells that don't appear to be replicating
    cn_g.to_csv(argv.output_g, index=False)

    # return the cells that look to be low quality
    cn_lq.to_csv(argv.output_lowqual, index=False)

