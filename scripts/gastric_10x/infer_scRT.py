from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.infer_scRT import scRT


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', help='input long-form copy number dataframe for all cells')
    p.add_argument('input_col', help='column that contains raw reads or rpm that is used as observed data in model')
    p.add_argument('clone_col', help='column in that contains phylogenetic clone id')
    p.add_argument('cn_col', help='column in that contains hmmcopy cn states')
    p.add_argument('gc_col', help='column containing gc values')
    p.add_argument('cn_prior_method', help='method for assigning the cn prior of each S-phase cell (i.e. g1_clones, g1_composite, diploid, etc)')
    p.add_argument('infer_mode', help='options: bulk/clone/cell/pyro')
    p.add_argument('cn_s_out', help='output tsv that is same as cn_s input with inferred scRT added')
    p.add_argument('supp_s_output', help='supplementerary output tsv containing sample- and library-level params inferred by the model')
    p.add_argument('cn_g_out', help='output tsv that is same as cn_g input with inferred scRT added')
    p.add_argument('supp_g_output', help='supplementerary output tsv containing sample- and library-level params inferred by the model')

    return p.parse_args()


def main():
    argv = get_args()

    # load data
    cn = pd.read_csv(argv.cn_input)

    # add a dummy column for library_id if it is not provided (i.e. dataset is one library)
    if 'library_id' not in cn.columns:
        cn['library_id'] = 1

    # split into initial G1/2 and S phase assignments based on the 'clonealign_clone_id' column
    cn_g = cn[cn[argv.clone_col]!='not_in_tree']
    cn_s = cn[cn[argv.clone_col]=='not_in_tree']


    # use library_id as the clone_id when it is not provided
    if argv.clone_col not in cn_g.columns:
        cn_g[argv.clone_col] = cn_g['library_id']

    # remove G1-phase clones containing <10 cells if using the composite prior
    if argv.infer_mode == 'g1_composite':
        counts = cn_g[['cell_id', argv.clone_col]].drop_duplicates().clone_id.value_counts()
        counts = counts.to_frame().reset_index()
        counts.columns = [argv.clone_col, 'num_cells']
        counts['freq'] = counts['num_cells'] / sum(counts['num_cells'])
        bad_clones = counts.query('num_cells < 10').clone_id.values
        cn_g = cn_g.loc[~cn_g[argv.clone_col].isin(bad_clones)]

    # temporarily remove columns that don't get used by infer_SPF in order to avoid
    # removing cells/loci that have NaN entries in some fields
    temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, 'library_id', argv.input_col]]
    temp_cn_g = cn_g[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, 'library_id', argv.clone_col, argv.input_col]]

    print('creating scrt object')
    # create SPF object with input
    scrt = scRT(temp_cn_s, temp_cn_g, input_col=argv.input_col, rt_prior_col=None, assign_col=argv.cn_col,
                cn_state_col=argv.cn_col, gc_col=argv.gc_col, cn_prior_method=argv.cn_prior_method, max_iter=1500)

    # run inference
    print('running inference')
    cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer_pyro_model()

    print('cn_s.shape', cn_s.shape)
    print('cn_s_with_scrt.shape', cn_s_with_scrt.shape)

    # merge cn_s_with_scrt with initial cn_s input to add columns that were excluded from temp_cn_s
    if argv.clone_col in cn_s_with_scrt.columns and argv.clone_col in cn_s.columns:
        cn_s_with_scrt.rename(columns={argv.clone_col: 'assigned_{}'.format(argv.clone_col)}, inplace=True)
    cn_s_out = pd.merge(cn_s, cn_s_with_scrt)
    print('cn_s_out.shape', cn_s_out.shape)
    cn_g_out = pd.merge(cn_g, cn_g_with_scrt)
    print('cn_g_out.shape', cn_g_out.shape)

    # save output files
    cn_s_out.to_csv(argv.cn_s_out, index=False)
    cn_g_out.to_csv(argv.cn_g_out, index=False)

    # this will be an empty df when argv.infer_mode!='pyro'
    supp_s_output.to_csv(argv.supp_s_output, index=False)
    supp_g_output.to_csv(argv.supp_g_output, index=False)


if __name__ == '__main__':
    main()
