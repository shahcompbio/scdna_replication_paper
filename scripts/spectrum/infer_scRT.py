from argparse import ArgumentParser
import pandas as pd
from scdna_replication_tools.infer_scRT import scRT


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells')
    p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells including clone_id')
    p.add_argument('rt_profile', help='input RT profile that should be used as a prior or initialization in the model')
    p.add_argument('input_col', help='column that contains raw reads or rpm that is used as observed data in model')
    p.add_argument('cn_col', help='column in that contains hmmcopy cn states')
    p.add_argument('copy_col', help='column in that contains hmmcopy cn states')
    p.add_argument('gc_col', help='column containing gc values')
    p.add_argument('cn_prior_method', help='method for assigning the cn prior of each S-phase cell (i.e. g1_clones, g1_composite, diploid, etc)')
    p.add_argument('infer_mode', help='options: bulk/clone/cell/pert')
    p.add_argument('rt_init_col', help='column containing the RT profile that should be used as initialization in the model')
    p.add_argument('cn_s_out', help='output csv that is same as cn_s input with inferred scRT added')
    p.add_argument('supp_s_output', help='supplementerary output csv containing sample- and library-level params inferred by the model')
    p.add_argument('cn_g_out', help='output csv that is same as cn_g input with inferred scRT added')
    p.add_argument('supp_g_output', help='supplementerary output csv containing sample- and library-level params inferred by the model')

    return p.parse_args()


def remove_hlamp_loci(cn_g, cn_s, argv, chrom=None, max_cn=11, pct_thresh=0.1):
    """
    Remove loci which contain high level amplifications `state==max_cn` and `copy>max_cn` in more than pct_thresh% of all the cells in cn_g.
    Restrict this filtering to just one chromosome if `chrom` is not None.
    """
    # compute the total number of cells in cn_g
    num_cells_g = len(cn_g['cell_id'].unique())
    # filter cn_s to all rows which have state==11 and copy>11
    hlamp_cn_g = cn_g.loc[(cn_g[argv.cn_col]==max_cn) & (cn_g[argv.copy_col]>max_cn)]
    # filter to just one chromosome if chrom is not None
    if chrom is not None:
        hlamp_cn_g = hlamp_cn_g.loc[hlamp_cn_g['chr']==chrom]
    # compute the number of cells each hlamp loci appears
    hlamp_cn_g = hlamp_cn_g.groupby(['chr', 'start', 'end']).cell_id.nunique().reset_index()
    # compute the fraction of cells each hlamp loci appears
    hlamp_cn_g['frac_cells'] = hlamp_cn_g['cell_id'] / num_cells_g
    # filter hlamp loci to those that appear in more than pct_thresh of all cells
    hlamp_cn_g = hlamp_cn_g.loc[hlamp_cn_g['frac_cells'] > pct_thresh]
    print('removing {} loci'.format(hlamp_cn_g.shape[0]))
    # remove all loci in cn_s and cn_g that appear in hlamp_cn_g
    merged_cn_g = pd.merge(cn_g, hlamp_cn_g[['chr', 'start', 'end']], how='outer', indicator=True)
    cn_g_filtered = merged_cn_g.loc[merged_cn_g['_merge']=='left_only'].drop(columns='_merge')
    merged_cn_s = pd.merge(cn_s, hlamp_cn_g[['chr', 'start', 'end']], how='outer', indicator=True)
    cn_s_filtered = merged_cn_s.loc[merged_cn_s['_merge']=='left_only'].drop(columns='_merge')
    return cn_g_filtered, cn_s_filtered


def main():
    argv = get_args()
    cn_s = pd.read_csv(argv.cn_s)
    cn_g = pd.read_csv(argv.cn_g1)
    rt_profile = pd.read_csv(argv.rt_profile)
    print('loaded data')

    # merge the rt profile into cn_s and cn_g
    print('before merging RT profile')
    print('cn_s.shape', cn_s.shape)
    print('cn_g.shape', cn_g.shape)
    print('rt_profile.shape', rt_profile.shape)
    cn_s = pd.merge(cn_s, rt_profile[['chr', 'start', 'end', argv.rt_init_col]])
    cn_g = pd.merge(cn_g, rt_profile[['chr', 'start', 'end', argv.rt_init_col]])
    print('after merging RT profile')
    print('cn_s.shape', cn_s.shape)
    print('cn_g.shape', cn_g.shape)

    # check if this is dataset SA1091
    dataset = argv.cn_s.split('/')[-2]
    if dataset == 'SA1091':
        print('removing hlamp loci')
        print('cn_g.shape', cn_g.shape)
        print('cn_s.shape', cn_s.shape)
        cn_g, cn_s = remove_hlamp_loci(cn_g, cn_s, argv)
        print('cn_g.shape', cn_g.shape) 
        print('cn_s.shape', cn_s.shape)
    

    # use library_id as the clone_id when it is not provided
    if 'clone_id' not in cn_g.columns:
        cn_g['clone_id'] = cn_g['library_id']

    # remove G1-phase clones containing <10 cells if using the composite prior
    if argv.cn_prior_method == 'g1_composite':
        counts = cn_g[['cell_id', 'clone_id']].drop_duplicates().clone_id.value_counts()
        counts = counts.to_frame().reset_index()
        counts.columns = ['clone_id', 'num_cells']
        counts['freq'] = counts['num_cells'] / sum(counts['num_cells'])
        bad_clones = counts.query('num_cells < 10').clone_id.values
        cn_g = cn_g.loc[~cn_g['clone_id'].isin(bad_clones)]

    # temporarily remove columns that don't get used by infer_SPF in order to avoid
    # removing cells/loci that have NaN entries in some fields
    temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.copy_col, 'library_id', argv.input_col, argv.rt_init_col]]
    temp_cn_g = cn_g[['cell_id', 'chr', 'start', 'end', argv.gc_col, argv.cn_col, argv.copy_col, 'library_id', 'clone_id', argv.input_col, argv.rt_init_col]]

    print('creating scrt object')
    # create scRT object with input
    scrt = scRT(
        temp_cn_s, temp_cn_g, input_col=argv.input_col, rt_init_col=argv.rt_init_col, assign_col=argv.copy_col,
        cn_state_col=argv.cn_col, gc_col=argv.gc_col, cn_prior_method=argv.cn_prior_method, max_iter=1500
    )

    # run inference
    print('running inference')
    cn_s_with_scrt, supp_s_output, cn_g_with_scrt, supp_g_output = scrt.infer(argv.infer_mode)

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
    cn_s_out.to_csv(argv.cn_s_out, index=False)
    cn_g_out.to_csv(argv.cn_g_out, index=False)

    # this will be an empty df when argv.infer_mode!='pyro'
    supp_s_output.to_csv(argv.supp_s_output, index=False)
    supp_g_output.to_csv(argv.supp_g_output, index=False)


if __name__ == '__main__':
    main()
