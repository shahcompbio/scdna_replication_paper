from argparse import ArgumentParser
import pandas as pd


def get_args():
	p = ArgumentParser()

	p.add_argument('--hmm', type=str, nargs='+', help='list of hmm reads data from different library and ticket ids')
	p.add_argument('--annotation', type=str, nargs='+', help='list of annotation metrics data from different library and ticket ids')
	p.add_argument('--samples', type=str, nargs='+', help='list of samples that should be included, all others (normal or controls) will be filtered out')
	p.add_argument('--output', type=str, help='output tsv file that combines all cells from same dataset into one file')

	return p.parse_args()


if __name__ == '__main__':
	argv = get_args()
	cn_pieces = []

	for f in argv.hmm:
		piece = pd.read_csv(
			f, index_col=['chr', 'start', 'end', 'width', 'cell_id'],
			# dtype=str
		)
		piece = piece[['reads', 'gc', 'map', 'copy', 'state']]
		cn_pieces.append(piece)

	cn = pd.concat(cn_pieces)
	cn = cn.reset_index()

	met_pieces = []
	for f in argv.annotation:
		piece = pd.read_csv(
			f, index_col=['cell_id']
		)
		piece = piece[['total_mapped_reads_hmmcopy', 'experimental_condition', 'is_s_phase', 'is_s_phase_prob',
					'quality', 'is_contaminated', 'cell_call', 'fastqscreen_grch37', 'fastqscreen_mm10', 'coverage_depth']]
		met_pieces.append(piece)

	metrics = pd.concat(met_pieces)
	print(metrics.columns)
	metrics = metrics.reset_index()
	metrics.rename(columns={'index': 'cell_id'}, inplace=True)

	metrics = metrics.loc[:, ~metrics.columns.duplicated()]

	# Remap experimental conditions and filter
	conditions = {
		'A': 'G1',
		'A-BSA': 'G1',
		'A-NCC': 'G1',
		'B': 'S',
		'B-NCC': 'S',
		'C': 'G2',
		'C-NCC': 'G2',
		'G1': 'G1',
		'G2': 'G2',
		'S': 'S',
		'D': 'D',
	}

	conditions = pd.Series(conditions)
	conditions.index.name = 'experimental_condition'
	conditions.name = 'cell_cycle_state'
	conditions = conditions.reset_index()
	metrics = metrics.drop_duplicates()
	print(conditions.head())
	print(metrics.columns)
	print(metrics.head())
	metrics = metrics.merge(conditions, how='left')
	metrics = metrics.drop_duplicates()

	cn = cn.merge(metrics, on='cell_id')
	cn = cn.loc[:, ~cn.columns.duplicated()]
	cn.drop_duplicates(inplace=True)

	# filter based on number of hmmcopy mapped reads
	cn = cn.query('total_mapped_reads_hmmcopy > 500000')

	# filter out control cells that don't contain one of the sample_ids in their cell_id
	print(argv.samples)
	cn = cn[cn['cell_id'].str.contains('|'.join(argv.samples))]

	# filter out cells based on quality, experimental condition and contamination	
	cn = cn[
		(cn['quality'] >= 0.75) &
		(cn['cell_call'].isin(['C1', 'C2']))
	]
	for rm_cond in ['gDNA', 'GM', 'NCC', 'NTC']:
		mask = ~cn['experimental_condition'].str.contains(rm_cond)
		cn = cn[mask]

	cn = cn[~cn['is_contaminated']]

	print("shape before filtering mouse bins", cn.shape)
	cn = cn[cn['fastqscreen_mm10'] < cn['fastqscreen_grch37']]
	print("shape after filtering mouse bins", cn.shape)

	# drop unmapable regions with negative gc values
	cn = cn[cn['gc']>=0]

	cn.to_csv(argv.output, sep='\t', index=False)
