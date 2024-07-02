from argparse import ArgumentParser
import pandas as pd
from scdna_replication_tools.infer_scRT import scRT


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells')
    p.add_argument('cn_g', help='input long-form copy number dataframe for G1/2-phase cells including clone_id')
    p.add_argument('input_col', help='column that contains raw reads or rpm that is used as observed data in model')
    p.add_argument('cn_col', help='column in that contains hmmcopy cn states')
    p.add_argument('copy_col', help='column in that contains hmmcopy cn states')
    p.add_argument('gc_col', help='column containing gc values')
    p.add_argument('cn_prior_method', help='method for assigning the cn prior of each S-phase cell (i.e. g1_clones, g1_composite, diploid, etc)')
    p.add_argument('cn_s_out', help='output tsv that is same as cn_s input with inferred scRT added')
    p.add_argument('supp_s_output', help='supplementerary output tsv containing sample- and library-level params inferred by the model')
    p.add_argument('cn_g_out', help='output tsv that is same as cn_g input with inferred scRT added')
    p.add_argument('supp_g_output', help='supplementerary output tsv containing sample- and library-level params inferred by the model')

    return p.parse_args()


def main():
    argv = get_args()
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g = pd.read_csv(argv.cn_g, sep='\t')
    print('loaded data')

    # use library_id as the clone_id when it is not provided
    if 'clone_id' not in cn_g.columns:
        cn_g['clone_id'] = cn_g['library_id']

    # temporarily remove columns that don't get used by infer_SPF in order to avoid
    # removing cells/loci that have NaN entries in some fields
    temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.copy_col, 'library_id', argv.input_col]]
    temp_cn_g = cn_g[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.copy_col, 'library_id', 'clone_id', argv.input_col]]

    print('creating scrt object')
    # create SPF object with input
    # run with composite cn prior to see what's going wrong
    scrt = scRT(temp_cn_s, temp_cn_g, input_col=argv.input_col, rt_prior_col=None, assign_col=argv.copy_col,
                cn_state_col=argv.cn_col, gc_col=argv.gc_col, cn_prior_method=argv.cn_prior_method, max_iter=1500)

    print('running inference')
    # run inference
    cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer_pyro_model()

    print('cn_s.shape', cn_s.shape)
    print('cn_s_with_scrt.shape', cn_s_with_scrt.shape)

    # merge cn_s_with_scrt with initial cn_s input to add columns that were excluded from temp_cn_s
    if 'clone_id' in cn_s_with_scrt.columns and 'clone_id' in cn_s.columns:
        cn_s_with_scrt.rename(columns={'clone_id': 'assigned_clone_id'}, inplace=True)
    cn_s_out = pd.merge(cn_s, cn_s_with_scrt)
    print('cn_s_out.shape', cn_s_out.shape)
    cn_g_out = pd.merge(cn_g, cn_g_with_scrt)
    print('cn_g_out.shape', cn_g_out.shape)

    # save output files
    cn_s_out.to_csv(argv.cn_s_out, sep='\t', index=False)
    cn_g_out.to_csv(argv.cn_g_out, sep='\t', index=False)

    # this will be an empty df when argv.infer_mode!='pyro'
    supp_s_output.to_csv(argv.supp_s_output, sep='\t', index=False)
    supp_g_output.to_csv(argv.supp_g_output, sep='\t', index=False)


if __name__ == '__main__':
    main()
