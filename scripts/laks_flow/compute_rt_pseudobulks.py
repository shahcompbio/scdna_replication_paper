import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('t47d_input', type=str, help='pyro model output for t47d s-phase cells')
    p.add_argument('GM18507_input', type=str, help='pyro model output for GM18507 s-phase cells')
    p.add_argument('merged_input', type=str, help='pyro model output for s-phase cells when run on both cell lines combined')
    p.add_argument('output', type=str, help='tsv file containing pseudobulk replication time for each locus')


    return p.parse_args()


def compute_loci_frac(cn):
    ''' Compute the fraction of replicated bins at each locus '''
    for (chrom, start), loci_cn in cn.groupby(['chr', 'start']):
        temp_rep = loci_cn['model_rep_state'].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[loci_cn.index, 'loci_frac_rep'] = temp_frac
    return cn


def compute_rt_pseudobulks(cn_t, cn_gm, cn_all):
    # compute the fraction of replicated bins at each locus
    cn_merged_t = compute_loci_frac(cn_all.query("sample_id=='SA1044'"))
    cn_merged_gm = compute_loci_frac(cn_all.query("sample_id=='SA928'"))
    cn_t = compute_loci_frac(cn_t)
    cn_gm = compute_loci_frac(cn_gm)
    print('computed loci frac')

    # create metric columns for each locus
    loci_metric_cols = ['chr', 'start', 'end', 'gc', 'loci_frac_rep']
    loci_metrics_merged_t = cn_merged_t[loci_metric_cols].drop_duplicates()
    loci_metrics_merged_gm = cn_merged_gm[loci_metric_cols].drop_duplicates()
    loci_metrics_t = cn_t[loci_metric_cols].drop_duplicates()
    loci_metrics_gm = cn_gm[loci_metric_cols].drop_duplicates()
    print('created metric columns for each locus')

    # rename columns based on merged vs split and cell line
    loci_metrics_merged_t = loci_metrics_merged_t.rename(columns={'loci_frac_rep': 'rt_merged_T47D'})
    loci_metrics_merged_gm = loci_metrics_merged_gm.rename(columns={'loci_frac_rep': 'rt_merged_GM18507'})
    loci_metrics_t = loci_metrics_t.rename(columns={'loci_frac_rep': 'rt_split_T47D'})
    loci_metrics_gm = loci_metrics_gm.rename(columns={'loci_frac_rep': 'rt_split_GM18507'})
    print('renamed columns')

    # merge together
    merged_loci_metrics = pd.merge(pd.merge(pd.merge(loci_metrics_merged_t, loci_metrics_merged_gm), loci_metrics_t), loci_metrics_gm)
    print('merged together')

    # compute the difference in RT between the two cell lines
    merged_loci_metrics['rt_diff_split'] = merged_loci_metrics['rt_split_T47D'] - merged_loci_metrics['rt_split_GM18507']
    merged_loci_metrics['rt_diff_merged'] = merged_loci_metrics['rt_merged_T47D'] - merged_loci_metrics['rt_merged_GM18507']

    print('returning merged copies')

    return merged_loci_metrics



if __name__=='__main__':
    argv = get_args()

    # load in data
    cn_t = pd.read_csv(argv.t47d_input, sep='\t')
    cn_gm = pd.read_csv(argv.GM18507_input, sep='\t')
    cn_all = pd.read_csv(argv.merged_input, sep='\t')

    print('data is loaded')

    bulk_rts = compute_rt_pseudobulks(cn_t, cn_gm, cn_all)

    bulk_rts.to_csv(argv.output, sep='\t', index=False)
    

