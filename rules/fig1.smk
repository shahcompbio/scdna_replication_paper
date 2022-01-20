import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"
samples = pd.read_csv('data/fitness/fitness_pseudobulk_qc_status.tsv', sep='\t')

# drop rows where ticket, library_id, or sample_id is NA
samples = samples[samples['mem_jira_ticket'].notna()]
samples = samples[samples['library_id'].notna()]
samples = samples[samples['sample_id'].notna()]

bad_datasets = []

rule all_fig1:
	input:
		expand(
			'analysis/fig1/{dataset}/spf_table.tsv',
			dataset=[
				d for d in config['fitness_datasets']
				if (d not in bad_datasets)
			]
		),
		expand(
			'plots/fig1/{dataset}/consensus_clone_copynumber.pdf',
			dataset=[
				d for d in config['fitness_datasets']
				if (d not in bad_datasets)
			]
		),
		expand(
			'plots/fig1/{dataset}/clonal_evolution.pdf',
			dataset=[
				d for d in config['fitness_datasets']
				if (d not in bad_datasets)
			]
		),
		expand(
			'plots/fig1/{dataset}/s_predictiveness.pdf',
			dataset=[
				d for d in config['fitness_datasets']
				if (d not in bad_datasets)
			]
		),
		


def dataset_cn_files(wildcards):
	mask = samples['datatag'] == wildcards.dataset
	library_ids = samples.loc[mask, 'library_id']
	ticket_ids = samples.loc[mask, 'mem_jira_ticket']

	files = expand(
		config['ticket_dir_500kb'] + \
			'/{ticket}/results/hmmcopy/{library}_reads.csv.gz',
		zip, ticket=ticket_ids, library=library_ids
	)
	return files


def dataset_sample_ids(wildcards):
	mask = samples['datatag'] == wildcards.dataset
	sample_ids = samples.loc[mask, 'sample_id']
	sample_ids = list(sample_ids)
	return sample_ids


def dataset_metric_files(wildcards):
	mask = samples['datatag'] == wildcards.dataset
	library_ids = samples.loc[mask, 'library_id']
	ticket_ids = samples.loc[mask, 'mem_jira_ticket']

	files = expand(
		config['ticket_dir_500kb'] + \
			'/{ticket}/results/annotation/{library}_metrics.csv.gz',
		zip, ticket=ticket_ids, library=library_ids
	)
	return files


rule collect_cn_data:
	input: 
		hmm = dataset_cn_files,
		annotation = dataset_metric_files
	output: 'analysis/fig1/{dataset}/{dataset}_cn_data.tsv'
	log: 'logs/fig1/{dataset}/collect_cn_data.log'
	params:
		samples = dataset_sample_ids
	shell: 
		'python scripts/fig1/collect_cn_data.py '
		'--hmm {input.hmm} --annotation {input.annotation} '
		'--samples {params.samples} --output {output} &> {log}'


rule assign_timepoints:
	input: 
		cn = 'analysis/fig1/{dataset}/{dataset}_cn_data.tsv',
		times = 'data/fitness/dlp_summaries_rebuttal.csv'
	output: 'analysis/fig1/{dataset}/{dataset}_cn_data_times.tsv'
	log: 'logs/fig1/{dataset}/assign_timepoints.log'
	shell:
		'python3 scripts/fig1/assign_timepoints.py '
		'{input} {output} &> {log}'


rule get_s_phase_cells:
	input: 'analysis/fig1/{dataset}/{dataset}_cn_data_times.tsv'
	output: 'analysis/fig1/{dataset}/s_phase_cells.tsv'
	run:
		df = pd.read_csv(str(input), sep='\t')
		df = df.query('is_s_phase_prob > 0.5')
		df.to_csv(str(output), sep='\t', index=False)


rule get_non_s_phase_cells:
	input:
		cn = 'analysis/fig1/{dataset}/{dataset}_cn_data_times.tsv',
		clones = 'data/fitness/fitness_cell_assignment_feb07_2020.tsv'
	output: 'analysis/fig1/{dataset}/non_s_phase_cells.tsv'
	run:
		df = pd.read_csv(str(input.cn), sep='\t')
		clones = pd.read_csv(str(input.clones), sep='\t')

		clones.rename(columns={'single_cell_id': 'cell_id',
							'letters': 'clone_id',
							'datatag': 'dataset_id'},
						inplace=True)

		# force clone cell_ids to match for the SA906 datasets
		clones['cell_id'] = clones['cell_id'].str.replace('SA906a', 'SA906', regex=False)
		clones['cell_id'] = clones['cell_id'].str.replace('SA906b', 'SA906', regex=False)

		df = df.query('is_s_phase_prob <= 0.5')
		# only use cells that have clone_id's assigned (and add clone_id column)
		df = pd.merge(df, clones)

		df.to_csv(str(output), sep='\t', index=False)


rule infer_SPF:
	input:
		cn_s = 'analysis/fig1/{dataset}/s_phase_cells.tsv',
		cn_g1 = 'analysis/fig1/{dataset}/non_s_phase_cells.tsv'
	output:
		cn_s_out = 'analysis/fig1/{dataset}/s_phase_cells_with_clones.tsv',
		spf_table = 'analysis/fig1/{dataset}/spf_table.tsv',
		clone_copy = 'analysis/fig1/{dataset}/clone_copy.tsv'
	params:
		input_col = 'copy'
	log: 'logs/fig1/{dataset}/infer_SPF.log'
	shell:
		'source ../scdna_replication_tools/venv/bin/activate ; '
		'python3 scripts/fig1/infer_SPF.py '
		'{input} {params} {output} &> {log} ; '
		'deactivate'


rule compute_consensus_clone_states:
	input: 'analysis/fig1/{dataset}/non_s_phase_cells.tsv'
	output: 'analysis/fig1/{dataset}/clone_states.tsv'
	params:
		input_col = 'state',
		clone_col = 'clone_id'
	log: 'logs/fig1/{dataset}/compute_consensus_clone_states.log'
	shell:
		'source ../scdna_replication_tools/venv/bin/activate ; '
		'python3 scripts/common/compute_consensus_clone_profiles.py '
		'{input} {params} {output} &> {log} ; '
		'deactivate'


rule plot_consensus_clone_copynumber:
	input:
		clone_states = 'analysis/fig1/{dataset}/clone_states.tsv',
		clone_copy = 'analysis/fig1/{dataset}/clone_copy.tsv'
	output: 'plots/fig1/{dataset}/consensus_clone_copynumber.pdf'
	log: 'logs/fig1/{dataset}/plot_consensus_clone_copynumber.log'
	shell:
		'source ../scgenome/venv/bin/activate ; '
		'python3 scripts/fig1/plot_consensus_clone_copynumber.py '
		'{input} {output} &> {log}'
		' ; deactivate'


rule plot_clonal_evolution:
	input: 
		s_phase = 'analysis/fig1/{dataset}/s_phase_cells_with_clones.tsv',
		non_s_phase = 'analysis/fig1/{dataset}/non_s_phase_cells.tsv'
	output:
		s_out = 'analysis/fig1/{dataset}/s_phase_clone_time_counts.tsv',
		non_s_out = 'analysis/fig1/{dataset}/non_s_phase_clone_time_counts.tsv',
		plot = 'plots/fig1/{dataset}/clonal_evolution.pdf'
	log: 'logs/fig1/{dataset}/plot_clonal_evolution.log'
	shell:
		'python3 scripts/fig1/plot_clonal_evolution.py '
		'{input} {output} &> {log}'


rule plot_s_predictiveness:
	input:
		s_phase = 'analysis/fig1/{dataset}/s_phase_clone_time_counts.tsv',
		non_s_phase = 'analysis/fig1/{dataset}/non_s_phase_clone_time_counts.tsv',
	output:
		tsv = 'analysis/fig1/{dataset}/s_predictiveness.tsv',
		plot = 'plots/fig1/{dataset}/s_predictiveness.pdf'
	params:
		dataset = lambda wildcards: wildcards.dataset
	log: 'logs/fig1/{dataset}/plot_s_predictiveness.log'
	shell:
		'python3 scripts/fig1/plot_s_predictiveness.py '
		'{input} {params} {output} &> {log}'