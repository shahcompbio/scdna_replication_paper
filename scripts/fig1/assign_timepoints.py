from argparse import ArgumentParser
import numpy as np
import pandas as pd

def get_args():
	p = ArgumentParser()

	p.add_argument('s_phase', help='input cell_id to clone_id tsv file for s phase cells')
	p.add_argument('non_s_phase', help='input cell_id to clone_id tsv file for non s phase cells')
	p.add_argument('time_legend', help='CSV file containing timepoints for each library')
	p.add_argument('s_out', help='output TSV that is same as s_phase input with timepoints added')
	p.add_argument('non_s_out', help='output TSV that is same as non_s_phase input with timepoints added')

	return p.parse_args()


def main():
	argv = get_args()
	s_df = pd.read_csv(argv.s_phase, sep='\t')
	non_s_df = pd.read_csv(argv.non_s_phase, sep='\t')
	time_df = pd.read_csv(argv.time_legend)

	# get rid of useless columns
	time_df = time_df[['library_id', 'label', 'datasetname', 'timepoint', 'sample_id']]

	# create library_id columns using cell_ids in the CN dataframes
	s_df['library_id'] = s_df['cell_id'].apply(lambda x: x.split('-')[1])
	non_s_df['library_id'] = non_s_df['cell_id'].apply(lambda x: x.split('-')[1])

	s_out = pd.merge(s_df, time_df)
	non_s_out = pd.merge(non_s_df, time_df)

	s_out.to_csv(argv.s_out, sep='\t')
	non_s_out.to_csv(argv.non_s_out, sep='\t')


if __name__ == '__main__':
	main()
