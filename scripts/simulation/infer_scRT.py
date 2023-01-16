from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.infer_scRT import scRT


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells')
    p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells including clone_id')
    p.add_argument('input_col', help='column in two cn dataframes to be used for matching S-phase cells to clones')
    p.add_argument('cn_col', help='column in that contains CN states for priors')
    p.add_argument('gc_col', help='column containing gc values')
    p.add_argument('cn_prior_method', help='method for assigning the cn prior of each S-phase cell (i.e. g1_clones, g1_composite, diploid, etc)')
    p.add_argument('infer_mode', help='options: bulk/clone/cell/pyro')
    p.add_argument('max_iter', type=int, help='max number of svi steps to take for each pyro model')
    p.add_argument('cn_s_out', help='output tsv that is same as cn_s input with inferred scRT added')
    p.add_argument('supp_s_output', help='supplementerary output tsv containing sample- and library-level params inferred by the model')
    p.add_argument('cn_g_out', help='output tsv that is same as cn_g input with inferred scRT added')
    p.add_argument('supp_g_output', help='supplementerary output tsv containing sample- and library-level params inferred by the model')

    return p.parse_args()


def main():
    argv = get_args()
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g1 = pd.read_csv(argv.cn_g1, sep='\t')
    print('loaded data')

    # make sure clone_id column is present if none are listed
    if 'clone_id' not in cn_g1.columns:
        cn_g1['clone_id'] = 'A'

    if 'library_id' not in cn_g1.columns:
        cn_g1['library_id'] = 'A'

    if 'library_id' not in cn_s.columns:
        cn_s['library_id'] = 'A'

    # temporarily remove columns that don't get used by infer_SPF in order to avoid
    # removing cells/loci that have NaN entries in some fields
    temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, 'library_id', argv.input_col]]
    temp_cn_g1 = cn_g1[['cell_id', 'chr', 'start', 'end', argv.gc_col, 'clone_id', argv.cn_col, 'library_id', argv.input_col]]

    print('creating scrt object')
    # create SPF object with input
    scrt = scRT(temp_cn_s, temp_cn_g1, input_col=argv.input_col, clone_col='clone_id', assign_col=argv.cn_col, rt_prior_col=None,
                cn_state_col=argv.cn_col, gc_col=argv.gc_col, cn_prior_method=argv.cn_prior_method, max_iter=argv.max_iter)

    # run inference
    print('running inference')
    cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer(level=argv.infer_mode)

    print('cn_s.shape', cn_s.shape)
    print('cn_s_with_scrt.shape', cn_s_with_scrt.shape)

    # merge cn_s_with_scrt with initial cn_s input to add columns that were excluded from temp_cn_s
    if 'clone_id' in cn_s_with_scrt.columns and 'clone_id' in cn_s.columns:
        cn_s_with_scrt.rename(columns={'clone_id': 'assigned_clone_id'}, inplace=True)
    cn_s_out = pd.merge(cn_s, cn_s_with_scrt)
    print('cn_s_out.shape', cn_s_out.shape)
    
    # only merge when cn_g_with_scrt is not empty
    if argv.infer_mode == 'pyro':
        cn_g_out = pd.merge(cn_g1, cn_g_with_scrt)
        print('cn_g_out.shape', cn_g_out.shape)
    else:
        cn_g_out = cn_g_with_scrt

    # save output files
    cn_s_out.to_csv(argv.cn_s_out, sep='\t', index=False)
    cn_g_out.to_csv(argv.cn_g_out, sep='\t', index=False)

    # this will be an empty df when argv.infer_mode!='pyro'
    supp_s_output.to_csv(argv.supp_s_output, sep='\t', index=False)
    supp_g_output.to_csv(argv.supp_g_output, sep='\t', index=False)


if __name__ == '__main__':
    main()
