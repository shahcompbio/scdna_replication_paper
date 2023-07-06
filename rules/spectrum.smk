import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = [
    'SPECTRUM-OV-002',
    'SPECTRUM-OV-014',
    'SPECTRUM-OV-068',
    'SPECTRUM-OV-105'
]

rule all_spectrum:
    input:
        expand(
            'plots/spectrum/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['spectrum_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/spectrum/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['spectrum_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/spectrum/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['spectrum_datasets']
                if (d not in bad_datasets)
            ]
        ),
         expand(
            'plots/spectrum/{dataset}/clone_spf.png',
            dataset=[
                d for d in config['spectrum_datasets']
                if (d not in bad_datasets)
            ]
        ),
         expand(
            'analysis/spectrum/{dataset}/scRT_pseudobulks.csv.gz',
            dataset=[
                d for d in config['spectrum_datasets']
                if (d not in bad_datasets)
            ]
        ),
         

rule collect_cn_data_sp:
    input: 'data/spectrum/experiment_analyses.csv'
    output: 'analysis/spectrum/{dataset}/cn_data.csv.gz'
    log: 'logs/spectrum/{dataset}/collect_cn_data.log'
    params:
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python3 scripts/spectrum/collect_cn_data.py '
        '--input {input} '
        '--dataset {params.dataset} '
        '--output {output} &> {log}'


rule clone_assignments_sp:
    input: 
        cn = 'analysis/spectrum/{dataset}/cn_data.csv.gz',
        clones = 'data/spectrum/clones/{dataset}_clones.csv'
    output: 'analysis/spectrum/{dataset}/cn_data_clones.csv.gz'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/spectrum/{dataset}/clone_assignments.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/clone_assignments.py '
        '{input} {params} {output} &> {log}'


rule compute_ccc_features_sp:
    input: 'analysis/spectrum/{dataset}/cn_data_clones.csv.gz'
    output: 'analysis/spectrum/{dataset}/cn_data_features.csv.gz'
    log: 'logs/spectrum/{dataset}/compute_ccc_features.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/compute_ccc_features.py '
        '{input} {params} {output} &> {log}'


rule plot_ccc_features_sp:
    input: 'analysis/spectrum/{dataset}/cn_data_features.csv.gz'
    output: 
        plot1 = 'plots/spectrum/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/spectrum/{dataset}/ccc_features_scatter.png'
    log: 'logs/spectrum/{dataset}/plot_ccc_features.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/plot_ccc_features.py '
        '{input} {params} {output} &> {log}'


rule split_cell_cycle_sp:
    input: 'analysis/spectrum/{dataset}/cn_data_features.csv.gz'
    output:
        cn_s = 'analysis/spectrum/{dataset}/s_phase_cells.csv.gz',
        cn_g1 = 'analysis/spectrum/{dataset}/g1_phase_cells.csv.gz'
    log: 'logs/spectrum/{dataset}/split_cell_cycle.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/split_cell_cycle.py '
        '{input} {params} {output} &> {log}'


rule plot_cn_heatmaps_sp:
    input:
        s_phase = 'analysis/spectrum/{dataset}/s_phase_cells.csv.gz',
        g1_phase = 'analysis/spectrum/{dataset}/g1_phase_cells.csv.gz'
    output: 'plots/spectrum/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/spectrum/{dataset}/plot_cn_heatmaps.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'


rule infer_scRT_pyro_sp:
    input:
        cn_s = 'analysis/spectrum/{dataset}/s_phase_cells.csv.gz',
        cn_g1 = 'analysis/spectrum/{dataset}/g1_phase_cells.csv.gz'
    output:
        main_s_out = 'analysis/spectrum/{dataset}/s_phase_cells_with_scRT.csv.gz',
        supp_s_out = 'analysis/spectrum/{dataset}/scRT_pyro_supp_s_output.csv.gz',
        main_g_out = 'analysis/spectrum/{dataset}/g1_phase_cells_with_scRT.csv.gz',
        supp_g_out = 'analysis/spectrum/{dataset}/scRT_pyro_supp_g_output.csv.gz'
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pert'
    log: 'logs/spectrum/{dataset}/infer_scRT_pyro.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/infer_scRT.py '
        '{input} {params} {output} &> {log}'


rule plot_pyro_model_output_sp:
    input:
        s_phase = 'analysis/spectrum/{dataset}/s_phase_cells_with_scRT.csv.gz',
        g1_phase = 'analysis/spectrum/{dataset}/g1_phase_cells_with_scRT.csv.gz'
    output: 'plots/spectrum/{dataset}/inferred_cn_rep_results.png',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        top_title_prefix = 'input_set:_unknown',
        bottom_title_prefix = 'input_set:_high-confidence_G1/2-phase'
    log: 'logs/spectrum/{dataset}/plot_pyro_model_output.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/plot_pyro_model_output.py '
        '--cn_s {input.s_phase} '
        '--cn_g {input.g1_phase} '
        '--dataset {params.dataset} '
        '--top_title_prefix {params.top_title_prefix} '
        '--bottom_title_prefix {params.bottom_title_prefix} '
        '--output {output} '
        '&> {log}'


rule predict_cycle_phase_sp:
    input: 
        cn_s = 'analysis/spectrum/{dataset}/s_phase_cells_with_scRT.csv.gz',
        cn_g = 'analysis/spectrum/{dataset}/g1_phase_cells_with_scRT.csv.gz'
    output: 
        out_s = 'analysis/spectrum/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        out_g = 'analysis/spectrum/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz',
        out_lowqual = 'analysis/spectrum/{dataset}/model_lowqual_cells.csv.gz'
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/spectrum/{dataset}/predict_cycle_phase.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/predict_cycle_phase.py '
        '{input} {params} {output} &> {log}'


rule plot_filtered_pyro_model_output_sp:
    input:
        s_phase = 'analysis/spectrum/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        g1_phase = 'analysis/spectrum/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz'
    output: 'plots/spectrum/{dataset}/inferred_cn_rep_results_filtered.png'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        top_title_prefix = 'PERT_S-phase',
        bottom_title_prefix = 'PERT_G1/2-phase'
    log: 'logs/spectrum/{dataset}/plot_filtered_pyro_model_output.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/plot_pyro_model_output.py '
        '--cn_s {input.s_phase} '
        '--cn_g {input.g1_phase} '
        '--dataset {params.dataset} '
        '--top_title_prefix {params.top_title_prefix} '
        '--bottom_title_prefix {params.bottom_title_prefix} '
        '--output {output} '
        '&> {log}'


rule compute_rt_pseudobulks_sp:
    input: 'analysis/spectrum/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz'
    output: 'analysis/spectrum/{dataset}/scRT_pseudobulks.csv.gz'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/spectrum/{dataset}/compute_rt_pseudobulks.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log}'


rule compute_cell_cycle_clone_counts_sp:
    input:
        cn_s = 'analysis/spectrum/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        cn_g = 'analysis/spectrum/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz'
    output: 'analysis/spectrum/{dataset}/cell_cycle_clone_counts.csv.gz'
    log: 'logs/spectrum/{dataset}/compute_cell_cycle_clone_counts.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/compute_cell_cycle_clone_counts.py '
        '{input} {output} &> {log}'


rule plot_clone_spf_sp:
    input: 'analysis/spectrum/{dataset}/cell_cycle_clone_counts.csv.gz'
    output: 'plots/spectrum/{dataset}/clone_spf.png'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/spectrum/{dataset}/plot_clone_spf.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/spectrum/plot_clone_spf.py '
        '{input} {params} {output} &> {log}'
