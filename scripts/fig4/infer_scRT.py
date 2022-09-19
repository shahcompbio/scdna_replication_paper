from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.infer_scRT import scRT


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells')
    p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells including clone_id')
    p.add_argument('cn_g2', help='input long-form copy number dataframe for G2-phase cells including clone_id')
    p.add_argument('input_col', help='column that contains raw reads or rpm that is used as observed data in model')
    p.add_argument('cn_col', help='column in that contains hmmcopy cn states')
    p.add_argument('copy_col', help='column in that contains hmmcopy cn states')
    p.add_argument('gc_col', help='column containing gc values')
    p.add_argument('cn_prior_method', help='method for assigning the cn prior of each S-phase cell (i.e. g1_clones, g1_composite, diploid, etc)')
    p.add_argument('infer_mode', help='options: bulk/clone/cell/pyro')
    p.add_argument('cn_s_out', help='output tsv that is same as cn_input with inferred scRT added')

    return p.parse_args()


def main():
    argv = get_args()
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g1 = pd.read_csv(argv.cn_g1, sep='\t')
    cn_g2 = pd.read_csv(argv.cn_g2, sep='\t')
    print('loaded data')

    # concate g1 and g2 phase cells together as they are both non-s phase
    cn_g = pd.concat([cn_g1, cn_g2], ignore_index=True)

    # use library_id as the clone_id when it is not provided
    if 'clone_id' not in cn_g.columns:
        cn_g['clone_id'] = cn_g['library_id']

    # temporarily remove columns that don't get used by infer_SPF in order to avoid
    # removing cells/loci that have NaN entries in some fields
    temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.copy_col, 'library_id', argv.input_col]]
    temp_cn_g = cn_g[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.copy_col, 'library_id', 'clone_id', argv.input_col]]

    print('creating scrt object')
    # create SPF object with input
    # running model with g1_cells prior method to see how model should work
    scrt = scRT(temp_cn_s, temp_cn_g, input_col=argv.input_col, rt_prior_col=None, assign_col=argv.copy_col,
                cn_state_col=argv.cn_col, gc_col=argv.gc_col, cn_prior_method='g1_clones')

    print('running inference')
    # run inference
    cn_s_with_scrt = scrt.infer_pyro_model(max_iter=20)
    print('done running inference')

    print('creating scrt object')
    # create SPF object with input
    # run with composite cn prior to see what's going wrong
    scrt = scRT(temp_cn_s, temp_cn_g, input_col=argv.input_col, rt_prior_col=None, assign_col=argv.copy_col,
                cn_state_col=argv.cn_col, gc_col=argv.gc_col, cn_prior_method=argv.cn_prior_method)

    print('running inference')
    # run inference
    cn_s_with_scrt = scrt.infer_pyro_model(max_iter=20)

    print('cn_s.shape', cn_s.shape)
    print('cn_s_with_scrt.shape', cn_s_with_scrt.shape)

    # merge cn_s_with_scrt with initial cn_s input to add columns that were excluded from temp_cn_s
    if 'clone_id' in cn_s_with_scrt.columns and 'clone_id' in cn_s.columns:
        cn_s_with_scrt.rename(columns={'clone_id': 'assigned_clone_id'}, inplace=True)
    cn_s_out = pd.merge(cn_s, cn_s_with_scrt)
    print('cn_s_out.shape', cn_s_out.shape)

    # save output files
    cn_s_out.to_csv(argv.cn_s_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
