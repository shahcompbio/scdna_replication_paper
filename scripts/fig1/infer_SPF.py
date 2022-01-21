from argparse import ArgumentParser
import numpy as np
import pandas as pd
from scdna_replication_tools.infer_SPF import SPF


def get_args():
	p = ArgumentParser()

	p.add_argument('cn_s', help='input long-form copy number dataframe for S-phase cells')
	p.add_argument('cn_g1', help='input long-form copy number dataframe for G1-phase cells including clone_id')
	p.add_argument('input_col', help='column in two cn dataframes to be used for matching S-phase cells to clones')
	p.add_argument('cn_s_out', help='output tsv that is same as cn_input with clone_id column added')
	p.add_argument('spf_table', help='table of S-phase fractions (with stdev) for each clone')
	p.add_argument('consensus_clones', help='consensus clone profiles used for assignment, matches input_col')

	return p.parse_args()


def main():
	argv = get_args()
	cn_s = pd.read_csv(argv.cn_s, sep='\t')
	cn_g1 = pd.read_csv(argv.cn_g1, sep='\t')
	print('loaded data')

	# temporarily remove columns that don't get used by infer_SPF in order to avoid
	# removing cells/loci that have NaN entries in some fields
	temp_cn_s = cn_s[['cell_id', 'chr', 'start', 'end', 'state', argv.input_col]]
	temp_cn_g1 = cn_g1[['cell_id', 'chr', 'start', 'end', 'clone_id', 'state', argv.input_col]]

	print('creating spf object')
	# create SPF object with input
	spf = SPF(temp_cn_s, temp_cn_g1, input_col=argv.input_col, clone_col='clone_id')

	print('running inference')
	# run inference
	cn_s_with_clone_id, spf_table = spf.infer()
	print('done running inference')

	# extract clone profiles from spf object
	clone_profiles = spf.clone_profiles

	print('cn_s.shape', cn_s.shape)
	print('cn_s_with_clone_id.shape', cn_s_with_clone_id.shape)

	# merge cn_s_with_clone_id with initial cn_s input to add columns that were excluded from temp_cn_s
	cn_s_out = pd.merge(cn_s, cn_s_with_clone_id)
	print('cn_s_out.shape', cn_s_out.shape)

	# save output files
	cn_s_out.to_csv(argv.cn_s_out, sep='\t', index=False)
	spf_table.to_csv(argv.spf_table, sep='\t', index=False)
	clone_profiles.to_csv(argv.consensus_clones, sep='\t')


if __name__ == '__main__':
	main()
