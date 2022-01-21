import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser


def get_args():
	p = ArgumentParser()

	p.add_argument('s_phase', help='input S phase cell counts for each time & clone')
	p.add_argument('non_s_phase', help='input non-S phase cell counts for each time & clone')
	p.add_argument('dataset')
	p.add_argument('output_tsv', help='output tsv file used to make the plot')
	p.add_argument('output_pdf', help='output pdf file fraction of each clone for S and non-S cells')

	return p.parse_args()


def main():
	argv = get_args()

	s_df = pd.read_csv(argv.s_phase, sep='\t')
	non_s_df = pd.read_csv(argv.non_s_phase, sep='\t')

	s_df['s_phase'] = True
	non_s_df['s_phase'] = False

	df = pd.concat([s_df, non_s_df]).reset_index().drop(columns=['index'])
	times = df.timepoint.unique()
	clones = df.clone_id.unique()
	df.set_index(['clone_id', 'timepoint', 's_phase'], inplace=True)

	time_diff_df = []
	phase_diff_df = []

	# find difference in a clone's number/fraction of cells between two adjacent timepoints
	for t in range(len(times)-1):
		for c in clones:
			t0 = times[t]
			t1 = times[t+1]
			for phase in [True, False]:
				count_diff = df.loc[c, t1, phase]['num_cells'] - df.loc[c, t0, phase]['num_cells']
				frac_diff = df.loc[c, t1, phase]['fraction'] - df.loc[c, t0, phase]['fraction']
				temp_time_df = pd.DataFrame({'clone_id': [c], 's_phase': [phase], 't0': [t0], 't1': [t1],
											'time_frac_diff': [frac_diff], 'time_count_diff': [count_diff]})
				time_diff_df.append(temp_time_df)

	time_diff_df = pd.concat(time_diff_df)

	# find difference between S-phase and G1-phase fractions for a given timepoint
	# large frac_diff values are situations with more S-phase cells than we'd expect
	for t in times:
		for c in clones:
			# how many cells are in S-phase for this time & clone compared to what we'd expect (non-S fraction)
			frac = df.loc[c, t, True]['fraction'] - df.loc[c, t, False]['fraction']
			temp_phase_df = pd.DataFrame({'clone_id': [c], 'timepoint': [t], 'phase_frac_diff': [frac]})
			phase_diff_df.append(temp_phase_df)

	phase_diff_df = pd.concat(phase_diff_df)
	
	phase_diff_df.set_index(['clone_id', 'timepoint'], inplace=True)
	out_df = time_diff_df.query('s_phase == False')
	out_df.drop(columns=['s_phase'], inplace=True)
	out_df['phase_frac_diff'] = None

	# for each time (t0) & clone, treat the phase_diff_df value as the "predicted" proliferation rate
	# and the time_diff_df value as the "observed" proliferation rate
	for i, row in out_df.iterrows():
		c = row.clone_id
		t = row.t0
		val = phase_diff_df.loc[c, t]['phase_frac_diff']
		out_df.loc[i, 'phase_frac_diff'] = val


	fig, ax = plt.subplots(1, 1, figsize=(6,6))

	sns.scatterplot(data=out_df, x='phase_frac_diff', y='time_frac_diff', hue='clone_id', ax=ax)
	ax.set(title=argv.dataset, xlabel='S minus G1 cell fraction of clone\nat time t0',
				ylabel='G1 cell fraction of clone\nat t1 minus t0')

	fig.savefig(argv.output_pdf, bbox_inches='tight')
	out_df.to_csv(argv.output_tsv, sep='\t', index=False)


if __name__ == '__main__':
	main()


