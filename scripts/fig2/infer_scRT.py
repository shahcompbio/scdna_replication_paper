from argparse import ArgumentParser
import numpy as np
import pandas as pd
import logging
from scdna_replication_tools.infer_scRT import scRT


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells')
    p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells including clone_id')
    p.add_argument('input_col', help='column in two cn dataframes to be used for matching S-phase cells to clones')
    p.add_argument('cn_s_out', help='output tsv that is same as cn_input with inferred scRT added')

    return p.parse_args()


def main():
    argv = get_args()
    cn_s = pd.read_csv(argv.cn_s, sep='\t')
    cn_g1 = pd.read_csv(argv.cn_g1, sep='\t')
    logging.info('loaded data')

    # make sure clone_id column is present if none are listed
    if 'clone_id' not in cn_g1.columns:
        cn_g1['clone_id'] = 'A'

    # use true G1 and rt states to recreate the 'state' column
    cn_g1['state'] = cn_g1['true_G1_state']
    cn_s['state'] = cn_s['true_G1_state'] * (1 + cn_s['true_rt_state'])

    # temporarily remove columns that don't get used by infer_SPF in order to avoid
    # removing cells/loci that have NaN entries in some fields
    temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', 'gc', 'state', argv.input_col]]
    temp_cn_g1 = cn_g1[['cell_id', 'chr', 'start', 'end', 'gc', 'clone_id', 'state', argv.input_col]]

    logging.info('creating scrt object')
    # create SPF object with input
    scrt = scRT(temp_cn_s, temp_cn_g1, input_col=argv.input_col, clone_col='clone_id')

    logging.info('running inference')
    # run inference
    cn_s_with_scrt = scrt.infer()
    logging.info('done running inference')

    logging.info('cn_s.shape {}'.format(cn_s.shape))
    logging.info('cn_s_with_clone_id.shape {}'.format(cn_s_with_clone_id.shape))

    # merge cn_s_with_clone_id with initial cn_s input to add columns that were excluded from temp_cn_s
    cn_s_out = pd.merge(cn_s, cn_s_with_scrt)
    logging.info('cn_s_out.shape {}'.format(cn_s_out.shape))

    # save output files
    cn_s_out.to_csv(argv.cn_s_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
