import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"
samples1 = pd.read_csv('data/fitness/dlp_summaries_rebuttal.csv')
samples2 = pd.read_csv('data/fitness/fitness_pseudobulk_qc_status.tsv', sep='\t')

# drop sample_id column from samples2
samples2.drop(columns=['sample_id'], inplace=True)

# merge samples1 and samples2 using 'library_id' to get samples
samples = pd.merge(samples1, samples2, on='library_id')

# drop rows where ticket, library_id, or sample_id is NA
samples = samples[samples['mem_jira_ticket'].notna()]
samples = samples[samples['library_id'].notna()]
samples = samples[samples['sample_id'].notna()]

bad_datasets = []

rule all_fig1_rx:
    input:
        expand(
            'plots/fig1_rx/{dataset}/cn_heatmaps.pdf',
            dataset=[
                d for d in config['fitness_rx_pairs']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig1_rx/{dataset}/clone_tree_heatmap.png',
            dataset=[
                d for d in config['fitness_rx_pairs']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig1_rx/{dataset}/consensus_clone_copynumber.pdf',
            dataset=[
                d for d in config['fitness_rx_pairs']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig1_rx/{dataset}/clonal_evolution_rx.pdf',
            dataset=[
                d for d in config['fitness_rx_pairs']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig1_rx/{dataset}/clonal_evolution_unrx.pdf',
            dataset=[
                d for d in config['fitness_rx_pairs']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig1_rx/{dataset}/total_SPF_vs_time.pdf',
            dataset=[
                d for d in config['fitness_rx_pairs']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig1_rx/{dataset}/first_unrx_SPF_vs_s_coeffs.pdf',
            dataset=[
                d for d in config['fitness_rx_pairs']
                if (d not in bad_datasets)
            ]
        ),


def dataset_cn_files(wildcards):
    rx_dataset = config['fitness_rx_pairs'][wildcards.dataset]['Rx']
    unrx_dataset = config['fitness_rx_pairs'][wildcards.dataset]['UnRx']
    rx_mask = samples['datasetname'] == rx_dataset
    unrx_mask = samples['datasetname'] == unrx_dataset
    
    library_ids = samples.loc[unrx_mask | rx_mask, 'library_id']
    ticket_ids = samples.loc[unrx_mask | rx_mask, 'mem_jira_ticket']

    files = expand(
        config['ticket_dir_500kb'] + \
            '/{ticket}/results/hmmcopy/{library}_reads.csv.gz',
        zip, ticket=ticket_ids, library=library_ids
    )

    return files


def dataset_sample_ids(wildcards):
    rx_dataset = config['fitness_rx_pairs'][wildcards.dataset]['Rx']
    unrx_dataset = config['fitness_rx_pairs'][wildcards.dataset]['UnRx']
    rx_mask = samples['datasetname'] == rx_dataset
    unrx_mask = samples['datasetname'] == unrx_dataset

    sample_ids = samples.loc[unrx_mask | rx_mask, 'sample_id']
    sample_ids = list(sample_ids)
    return sample_ids


def dataset_metric_files(wildcards):
    rx_dataset = config['fitness_rx_pairs'][wildcards.dataset]['Rx']
    unrx_dataset = config['fitness_rx_pairs'][wildcards.dataset]['UnRx']
    rx_mask = samples['datasetname'] == rx_dataset
    unrx_mask = samples['datasetname'] == unrx_dataset
    
    library_ids = samples.loc[unrx_mask | rx_mask, 'library_id']
    ticket_ids = samples.loc[unrx_mask | rx_mask, 'mem_jira_ticket']

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
    output: 'analysis/fig1_rx/{dataset}/{dataset}_cn_data.tsv'
    log: 'logs/fig1_rx/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids
    shell: 
        'python scripts/fig1_rx/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --output {output} &> {log}'


rule assign_timepoints:
    input: 
        cn = 'analysis/fig1_rx/{dataset}/{dataset}_cn_data.tsv',
        times = 'data/fitness/dlp_summaries_rebuttal.csv'
    output: 'analysis/fig1_rx/{dataset}/{dataset}_cn_data_times.tsv'
    log: 'logs/fig1_rx/{dataset}/assign_timepoints.log'
    shell:
        'python3 scripts/fig1_rx/assign_timepoints.py '
        '{input} {output} &> {log}'


rule get_s_phase_cells:
    input: 'analysis/fig1_rx/{dataset}/{dataset}_cn_data_times.tsv'
    output: 'analysis/fig1_rx/{dataset}/s_phase_cells.tsv'
    run:
        df = pd.read_csv(str(input), sep='\t')
        df = df.query('is_s_phase_prob > 0.5')
        df.to_csv(str(output), sep='\t', index=False)


rule get_non_s_phase_cells:
    input:
        cn = 'analysis/fig1_rx/{dataset}/{dataset}_cn_data_times.tsv',
        clones = 'data/fitness/fitness_cell_assignment_feb07_2020.tsv'
    output: 'analysis/fig1_rx/{dataset}/non_s_phase_cells.tsv'
    run:
        df = pd.read_csv(str(input.cn), sep='\t')
        clones = pd.read_csv(str(input.clones))

        clones.rename(columns={'single_cell_id': 'cell_id',
                            'letters': 'clone_id',
                            'datatag': 'dataset_id'},
                        inplace=True)
        clones.drop(columns=['sample_id', 'V1'], inplace=True)

        # force clone cell_ids to match for the SA906 datasets
        clones['cell_id'] = clones['cell_id'].str.replace('SA906a', 'SA906', regex=False)
        clones['cell_id'] = clones['cell_id'].str.replace('SA906b', 'SA906', regex=False)

        df = df.query('is_s_phase_prob <= 0.5')
        # only use cells that have clone_id's assigned (and add clone_id column)
        df = pd.merge(df, clones)

        df.to_csv(str(output), sep='\t', index=False)


rule infer_SPF:
    input:
        cn_s = 'analysis/fig1_rx/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig1_rx/{dataset}/non_s_phase_cells.tsv'
    output:
        cn_s_out = 'analysis/fig1_rx/{dataset}/s_phase_cells_with_clones.tsv',
        spf_table = 'analysis/fig1_rx/{dataset}/spf_table.tsv',
        clone_copy = 'analysis/fig1_rx/{dataset}/clone_copy.tsv'
    params:
        input_col = 'copy'
    log: 'logs/fig1_rx/{dataset}/infer_SPF.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig1_rx/infer_SPF.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_consensus_clone_states:
    input: 'analysis/fig1_rx/{dataset}/non_s_phase_cells.tsv'
    output: 'analysis/fig1_rx/{dataset}/clone_states.tsv'
    params:
        input_col = 'state',
        clone_col = 'clone_id'
    log: 'logs/fig1_rx/{dataset}/compute_consensus_clone_states.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/common/compute_consensus_clone_profiles.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_consensus_clone_copynumber:
    input:
        clone_states = 'analysis/fig1_rx/{dataset}/clone_states.tsv',
        clone_copy = 'analysis/fig1_rx/{dataset}/clone_copy.tsv'
    output: 'plots/fig1_rx/{dataset}/consensus_clone_copynumber.pdf'
    log: 'logs/fig1_rx/{dataset}/plot_consensus_clone_copynumber.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig1_rx/plot_consensus_clone_copynumber.py '
        '{input} {output} &> {log}'
        ' ; deactivate'


rule plot_clone_tree_heatmap:
    input: 'analysis/fig1_rx/{dataset}/clone_states.tsv',
    output: 'plots/fig1_rx/{dataset}/clone_tree_heatmap.png'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig1_rx/{dataset}/plot_clone_tree_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig1_rx/plot_clone_tree_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_cn_heatmaps:
    input:
        s_phase = 'analysis/fig1_rx/{dataset}/s_phase_cells_with_clones.tsv',
        non_s_phase = 'analysis/fig1_rx/{dataset}/non_s_phase_cells.tsv'
    output: 'plots/fig1_rx/{dataset}/cn_heatmaps.pdf'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig1_rx/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig1_rx/plot_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'

rule split_by_rx:
    input:
        s_phase = 'analysis/fig1_rx/{dataset}/s_phase_cells_with_clones.tsv',
        non_s_phase = 'analysis/fig1_rx/{dataset}/non_s_phase_cells.tsv'
    output:
        s_phase_rx = 'analysis/fig1_rx/{dataset}/s_phase_rx_cells_with_clones.tsv',
        non_s_phase_rx = 'analysis/fig1_rx/{dataset}/non_s_phase_rx_cells.tsv',
        s_phase_unrx = 'analysis/fig1_rx/{dataset}/s_phase_unrx_cells_with_clones.tsv',
        non_s_phase_unrx = 'analysis/fig1_rx/{dataset}/non_s_phase_unrx_cells.tsv'
    params:
        rx_datasetname = lambda wildcards: config['fitness_rx_pairs'][wildcards.dataset]['Rx'],
        unrx_datasetname = lambda wildcards: config['fitness_rx_pairs'][wildcards.dataset]['UnRx']
    run:
        s_df = pd.read_csv(str(input.s_phase), sep='\t')
        g_df = pd.read_csv(str(input.non_s_phase), sep='\t')

        s_rx_df = s_df.loc[s_df['datasetname']==str(params.rx_datasetname)]
        g_rx_df = g_df.loc[g_df['datasetname']==str(params.rx_datasetname)]
        s_unrx_df = s_df.loc[s_df['datasetname']==str(params.unrx_datasetname)]
        g_unrx_df = g_df.loc[g_df['datasetname']==str(params.unrx_datasetname)]

        s_rx_df.to_csv(str(output.s_phase_rx), sep='\t', index=False)
        g_rx_df.to_csv(str(output.non_s_phase_rx), sep='\t', index=False)
        s_unrx_df.to_csv(str(output.s_phase_unrx), sep='\t', index=False)
        g_unrx_df.to_csv(str(output.non_s_phase_unrx), sep='\t', index=False)


rule plot_clonal_evolution_unrx:
    input: 
        s_phase = 'analysis/fig1_rx/{dataset}/s_phase_unrx_cells_with_clones.tsv',
        non_s_phase = 'analysis/fig1_rx/{dataset}/non_s_phase_unrx_cells.tsv',
        times = 'data/fitness/fitness_time_scale.tsv'
    output:
        s_out = 'analysis/fig1_rx/{dataset}/s_phase_unrx_clone_time_counts.tsv',
        non_s_out = 'analysis/fig1_rx/{dataset}/non_s_phase_unrx_clone_time_counts.tsv',
        plot = 'plots/fig1_rx/{dataset}/clonal_evolution_unrx.pdf'
    params:
        dataset = lambda wildcards: config['fitness_rx_pairs'][wildcards.dataset]['UnRx']
    log: 'logs/fig1_rx/{dataset}/plot_clonal_evolution.log'
    shell:
        'python3 scripts/fig1_rx/plot_clonal_evolution.py '
        '{input} {params} {output} &> {log}'


rule plot_clonal_evolution_rx:
    input: 
        s_phase = 'analysis/fig1_rx/{dataset}/s_phase_rx_cells_with_clones.tsv',
        non_s_phase = 'analysis/fig1_rx/{dataset}/non_s_phase_rx_cells.tsv',
        times = 'data/fitness/fitness_time_scale.tsv'
    output:
        s_out = 'analysis/fig1_rx/{dataset}/s_phase_rx_clone_time_counts.tsv',
        non_s_out = 'analysis/fig1_rx/{dataset}/non_s_phase_rx_clone_time_counts.tsv',
        plot = 'plots/fig1_rx/{dataset}/clonal_evolution_rx.pdf'
    params:
        dataset = lambda wildcards: config['fitness_rx_pairs'][wildcards.dataset]['Rx']
    log: 'logs/fig1_rx/{dataset}/plot_clonal_evolution.log'
    shell:
        'python3 scripts/fig1_rx/plot_clonal_evolution.py '
        '{input} {params} {output} &> {log}'


rule plot_total_SPF_vs_time:
    input:
        s_phase_rx = 'analysis/fig1_rx/{dataset}/s_phase_rx_clone_time_counts.tsv',
        non_s_phase_rx = 'analysis/fig1_rx/{dataset}/non_s_phase_rx_clone_time_counts.tsv',
        s_phase_unrx = 'analysis/fig1_rx/{dataset}/s_phase_unrx_clone_time_counts.tsv',
        non_s_phase_unrx = 'analysis/fig1_rx/{dataset}/non_s_phase_unrx_clone_time_counts.tsv'
    output: 'plots/fig1_rx/{dataset}/total_SPF_vs_time.pdf'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig1_rx/{dataset}/plot_total_SPF_vs_time.log'
    shell:
        'python3 scripts/fig1_rx/plot_total_SPF_vs_time.py '
        '{input} {params} {output} &> {log}'


rule plot_first_unrx_SPF_vs_s_coeffs:
    input:
        s_phase_rx = 'analysis/fig1_rx/{dataset}/s_phase_rx_clone_time_counts.tsv',
        non_s_phase_rx = 'analysis/fig1_rx/{dataset}/non_s_phase_rx_clone_time_counts.tsv',
        s_phase_unrx = 'analysis/fig1_rx/{dataset}/s_phase_unrx_clone_time_counts.tsv',
        non_s_phase_unrx = 'analysis/fig1_rx/{dataset}/non_s_phase_unrx_clone_time_counts.tsv',
        s_coeffs = 'data/fitness/fitness_selection_coefficients.csv'
    output: 
        pdf = 'plots/fig1_rx/{dataset}/first_unrx_SPF_vs_s_coeffs.pdf',
        tsv = 'analysis/fig1_rx/{dataset}/first_unrx_SPF_vs_s_coeffs.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        rx_name = lambda wildcards: config['fitness_rx_pairs'][wildcards.dataset]['Rx'],
        unrx_name = lambda wildcards: config['fitness_rx_pairs'][wildcards.dataset]['UnRx']
    log: 'logs/fig1_rx/{dataset}/plot_first_unrx_SPF_vs_s_coeffs.log'
    shell:
        'python3 scripts/fig1_rx/plot_first_unrx_SPF_vs_s_coeffs.py '
        '{input} {params} {output} &> {log}'


