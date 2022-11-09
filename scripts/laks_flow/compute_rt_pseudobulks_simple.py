from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles


def get_args():
    p = ArgumentParser()

    p.add_argument('input', help='long-form copy number dataframe for S-phase cells with scRT data')
    p.add_argument('rep_col', type=str, help='column name for replicated status of each bin in pyro model')
    p.add_argument('output', help='scRT pseudobulks')

    return p.parse_args()


def compute_loci_frac(cn, argv):
    ''' Compute the fraction of replicated bins at each locus '''
    for (chrom, start), loci_cn in cn.groupby(['chr', 'start']):
        temp_rep = loci_cn[argv.rep_col].values
        temp_frac = sum(temp_rep) / len(temp_rep)
        cn.loc[loci_cn.index, 'loci_frac_rep'] = temp_frac
    return cn


def compute_rt_pseudobulks(cn, argv):
    # compute the fraction of replicated bins at each locus
    cn_t = compute_loci_frac(cn.query("sample_id=='SA1044'"), argv)
    cn_gm = compute_loci_frac(cn.query("sample_id=='SA928'"), argv)
    cn_all = compute_loci_frac(cn, argv)

    # create metric columns for each locus
    loci_metric_cols = ['chr', 'start', 'end', 'gc', 'loci_frac_rep']
    loci_metrics_t = cn_t[loci_metric_cols].drop_duplicates()
    loci_metrics_gm = cn_gm[loci_metric_cols].drop_duplicates()
    loci_metrics_all = cn_all[loci_metric_cols].drop_duplicates()

    # rename columns based on joint vs split and cell line
    loci_metrics_t = loci_metrics_t.rename(columns={'loci_frac_rep': 'rt_T47D'})
    loci_metrics_gm = loci_metrics_gm.rename(columns={'loci_frac_rep': 'rt_GM18507'})
    loci_metrics_all = loci_metrics_all.rename(columns={'loci_frac_rep': 'rt_all'})

    # merge together
    merged_loci_metrics = pd.merge(pd.merge(loci_metrics_t, loci_metrics_gm), loci_metrics_all)

    # compute the difference in RT between the two cell lines
    merged_loci_metrics['rt_diff'] = merged_loci_metrics['rt_T47D'] - merged_loci_metrics['rt_GM18507']

    return merged_loci_metrics


def main():
    argv = get_args()

    df = pd.read_csv(argv.input, sep='\t')

    bulk_rts = compute_rt_pseudobulks(df, argv)

    bulk_rts.to_csv(argv.output, sep='\t', index=False)


if __name__ == '__main__':
    main()