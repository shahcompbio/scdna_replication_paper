import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"
hmmcopy_samples = pd.read_csv('data/signatures/signatures-hmmcopy.csv')
metrics_samples = pd.read_csv('data/signatures/signatures-annotation.csv')
htert_hmmcopy = hmmcopy_samples.loc[hmmcopy_samples['isabl_patient_id']=='184-hTERT']
htert_metrics = metrics_samples.loc[metrics_samples['isabl_patient_id']=='184-hTERT']

bad_datasets = []

rule all_fitness_lines:
    input:
        expand(
            'plots/fitness_lines/{dataset}/ccc_features_scatter.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/clone_spf.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/clone_rt.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/cn_pseudobulks1.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/frac_rt_distributions.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/clones_vs_time.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        'plots/fitness_lines/s_predictiveness.png'


def dataset_cn_files_fl(wildcards):
    if wildcards.dataset == 'OV2295':
        files = hmmcopy_samples.loc[
            hmmcopy_samples['isabl_patient_id']==wildcards.dataset].loc[
            hmmcopy_samples['result_type']=='reads']['result_filepath'].values
    else:
        col = 'in_{}'.format(wildcards.dataset)
        htert_hmmcopy[col] = list(
            map(lambda x: x.startswith(wildcards.dataset), htert_hmmcopy['isabl_sample_id'])
        )
        files = htert_hmmcopy.loc[
            htert_hmmcopy[col]==True].loc[
            htert_hmmcopy['result_type']=='reads']['result_filepath'].values
    
    return expand(files)


def dataset_metric_files_fl(wildcards):
    if wildcards.dataset == 'OV2295':
        files = metrics_samples.loc[
            metrics_samples['isabl_patient_id']==wildcards.dataset].loc[
            metrics_samples['result_type']=='metrics']['result_filepath'].values
    else:
        col = 'in_{}'.format(wildcards.dataset)
        htert_metrics[col] = list(
            map(lambda x: x.startswith(wildcards.dataset), htert_metrics['isabl_sample_id'])
        )
        files = htert_metrics.loc[
            htert_metrics[col]==True].loc[
            htert_metrics['result_type']=='metrics']['result_filepath'].values
    
    return expand(files)


def dataset_sample_ids_fl(wildcards):
    if wildcards.dataset == 'OV2295':
        sample_ids = metrics_samples.loc[metrics_samples['isabl_patient_id']==wildcards.dataset]['isabl_sample_id'].unique()
    else:
        col = 'in_{}'.format(wildcards.dataset)
        htert_metrics[col] = list(
            map(lambda x: x.startswith(wildcards.dataset), htert_metrics['isabl_sample_id'])
        )
        sample_ids = htert_metrics.loc[htert_metrics[col]==True]['isabl_sample_id'].unique()
    
    return expand(sample_ids)


rule collect_cn_data_fl:
    input: 
        hmm = dataset_cn_files_fl,
        annotation = dataset_metric_files_fl
    output: 'analysis/fitness_lines/{dataset}/cn_data.tsv'
    log: 'logs/fitness_lines/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids_fl,
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python scripts/fitness_lines/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --dataset {params.dataset} '
        '--output {output} &> {log}'


rule clone_assignments_fl:
    input: 
        cn = 'analysis/fitness_lines/{dataset}/cn_data.tsv',
        clones = 'data/fitness/fitness_cell_assignment_feb07_2020.tsv'
    output: 'analysis/fitness_lines/{dataset}/cn_data_clones.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/fitness_lines/{dataset}/clone_assignments.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/clone_assignments.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_ccc_features_fl:
    input: 'analysis/fitness_lines/{dataset}/cn_data_clones.tsv'
    output: 'analysis/fitness_lines/{dataset}/cn_data_features.tsv'
    log: 'logs/fitness_lines/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_fl:
    input: 'analysis/fitness_lines/{dataset}/cn_data_features.tsv'
    output: 
        plot1 = 'plots/fitness_lines/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/fitness_lines/{dataset}/ccc_features_scatter.png'
    log: 'logs/fitness_lines/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule split_cell_cycle_fl:
    input: 'analysis/fitness_lines/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fitness_lines/{dataset}/g1_phase_cells.tsv'
    log: 'logs/fitness_lines/{dataset}/split_cell_cycle.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/split_cell_cycle.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_fl:
    input:
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fitness_lines/{dataset}/g1_phase_cells.tsv'
    output:
        main_s_out = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        supp_s_out = 'analysis/fitness_lines/{dataset}/scRT_pyro_supp_s_output.tsv',
        main_g_out = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT.tsv',
        supp_g_out = 'analysis/fitness_lines/{dataset}/scRT_pyro_supp_g_output.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pyro'
    log: 'logs/fitness_lines/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_pyro_model_output_fl:
    input:
        s_phase = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        g1_phase = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT.tsv'
    output:
        plot1 = 'plots/fitness_lines/{dataset}/inferred_cn_rep_results.png',
        plot2 = 'plots/fitness_lines/{dataset}/s_vs_g_hmmcopy_states.png',
        plot3 = 'plots/fitness_lines/{dataset}/s_vs_g_rpm.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule revise_cell_cycle_labels_fl:
    input: 
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        cn_g = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT.tsv',
    output:
        out_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered_no_times.tsv',
        out_g = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT_filtered_no_times.tsv',
        out_lowqual = 'analysis/fitness_lines/{dataset}/model_lowqual_cells_no_times.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/fitness_lines/{dataset}/revise_cell_cycle_labels.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/revise_cell_cycle_labels.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# use the model output files from sig_lines.smk to assign S-phase cells to timepoints
rule assign_timepoints_fl:
    input:
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered_no_times.tsv',
        cn_g = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT_filtered_no_times.tsv',
        cn_lowqual = 'analysis/fitness_lines/{dataset}/model_lowqual_cells_no_times.tsv',
        times = 'data/fitness/dlp_summaries_rebuttal.csv'
    output: 
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        cn_g = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        cn_lowqual = 'analysis/fitness_lines/{dataset}/model_lowqual_cells.tsv',
    log: 'logs/fitness_lines/{dataset}/assign_timepoints.log'
    shell:
        'python3 scripts/fitness_lines/assign_timepoints.py '
        '{input} {output} &> {log}'


rule plot_filtered_pyro_model_output_fl:
    input:
        s_phase = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv'
    output:
        plot1 = 'plots/fitness_lines/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/fitness_lines/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/fitness_lines/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/common/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_filtered_fl:
    input: 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 
        plot1 = 'plots/fitness_lines/{dataset}/cn_scRT_heatmaps.png',
        plot2 = 'plots/fitness_lines/{dataset}/frac_rt_distributions.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/fitness_lines/{dataset}/plot_inferred_cn_vs_scRT_filtered.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_fl:
    input: 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/fitness_lines/{dataset}/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/fitness_lines/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_cn_pseudobulks_fl:
    input: 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/fitness_lines/{dataset}/cn_pseudobulks.tsv'
    params:
        cn_state_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/compute_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/compute_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_fl:
    input:
        s_phase = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv'
    output: 'plots/fitness_lines/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fitness_lines/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_clone_rt_and_spf_fl:
    input: 
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        cn_g = 'analysis/fitness_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        rt = 'analysis/fitness_lines/{dataset}/scRT_pseudobulks.tsv'
    output:
        tsv = 'analysis/fitness_lines/{dataset}/cell_cycle_clone_counts.tsv',
        clone_rt = 'plots/fitness_lines/{dataset}/clone_rt.png',
        clone_spf = 'plots/fitness_lines/{dataset}/clone_spf.png'
    params:
        rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_clone_rt_and_spf.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/plot_clone_rt_and_spf.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_clones_vs_time_fl:
    input: 'analysis/fitness_lines/{dataset}/cell_cycle_clone_counts.tsv'
    output: 
        plot1 = 'plots/fitness_lines/{dataset}/clones_vs_time.png',
        plot2 = 'plots/fitness_lines/{dataset}/total_cells_vs_time.png'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_clones_vs_time.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/plot_clones_vs_time.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_pseudobulks_fl:
    input: 'analysis/fitness_lines/{dataset}/cn_pseudobulks.tsv'
    output: 
        plot1 = 'plots/fitness_lines/{dataset}/cn_pseudobulks1.png',
        plot2 = 'plots/fitness_lines/{dataset}/cn_pseudobulks2.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/plot_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule s_predictiveness_fl:
    input:
        SA039 = 'analysis/fitness_lines/SA039/cell_cycle_clone_counts.tsv',
        SA906a = 'analysis/fitness_lines/SA906a/cell_cycle_clone_counts.tsv',
        SA906b = 'analysis/fitness_lines/SA906b/cell_cycle_clone_counts.tsv',
    output:
        tsv = 'analysis/fitness_lines/s_predictiveness.tsv',
        plot = 'plots/fitness_lines/s_predictiveness.png'
    log: 'logs/fitness_lines/s_predictiveness.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/fitness_lines/s_predictiveness.py '
        '{input} {output} &> {log} ; '
        'deactivate'