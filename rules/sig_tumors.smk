import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = ['SA535', 'SA609']

rule all_sig_tumors:
    input:
        expand(
            'plots/sig_tumors/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/twidth_curves.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/clone_spf.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/clone_rt.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        'plots/sig_tumors/cohort_spf.png'
        


rule collect_cn_data_st:
    input: 
        hmm = 'data/signatures/signatures-hmmcopy.csv',
        annotation = 'data/signatures/signatures-annotation.csv'
    output: 'analysis/sig_tumors/{dataset}/cn_data.csv.gz'
    log: 'logs/sig_tumors/{dataset}/collect_cn_data.log'
    params:
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python scripts/sig_tumors/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--dataset {params.dataset} '
        '--output {output} &> {log}'


rule clone_assignments_st:
    input: 
        cn = 'analysis/sig_tumors/{dataset}/cn_data.csv.gz',
        clones = 'data/signatures/clone_trees/{dataset}_clones.tsv'
    output: 'analysis/sig_tumors/{dataset}/cn_data_clones.csv.gz'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/sig_tumors/{dataset}/clone_assignments.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/clone_assignments.py '
        '{input} {params} {output} &> {log}'


rule compute_ccc_features_st:
    input: 'analysis/sig_tumors/{dataset}/cn_data_clones.csv.gz'
    output: 'analysis/sig_tumors/{dataset}/cn_data_features.csv.gz'
    log: 'logs/sig_tumors/{dataset}/compute_ccc_features.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/compute_ccc_features.py '
        '{input} {params} {output} &> {log}'


rule plot_ccc_features_st:
    input: 'analysis/sig_tumors/{dataset}/cn_data_features.csv.gz'
    output: 
        plot1 = 'plots/sig_tumors/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/sig_tumors/{dataset}/ccc_features_scatter.png'
    log: 'logs/sig_tumors/{dataset}/plot_ccc_features.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/plot_ccc_features.py '
        '{input} {params} {output} &> {log}'


rule split_cell_cycle_st:
    input: 'analysis/sig_tumors/{dataset}/cn_data_features.csv.gz'
    output:
        cn_s = 'analysis/sig_tumors/{dataset}/s_phase_cells.csv.gz',
        cn_g1 = 'analysis/sig_tumors/{dataset}/g1_phase_cells.csv.gz'
    log: 'logs/sig_tumors/{dataset}/split_cell_cycle.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/split_cell_cycle.py '
        '{input} {params} {output} &> {log}'


rule plot_cn_heatmaps_st:
    input:
        s_phase = 'analysis/sig_tumors/{dataset}/s_phase_cells.csv.gz',
        g1_phase = 'analysis/sig_tumors/{dataset}/g1_phase_cells.csv.gz'
    output: 'plots/sig_tumors/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_cn_heatmaps.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'


rule infer_scRT_pyro_st:
    input:
        cn_s = 'analysis/sig_tumors/{dataset}/s_phase_cells.csv.gz',
        cn_g1 = 'analysis/sig_tumors/{dataset}/g1_phase_cells.csv.gz'
    output:
        main_s_out = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT.csv.gz',
        supp_s_out = 'analysis/sig_tumors/{dataset}/scRT_pyro_supp_s_output.csv.gz',
        main_g_out = 'analysis/sig_tumors/{dataset}/g1_phase_cells_with_scRT.csv.gz',
        supp_g_out = 'analysis/sig_tumors/{dataset}/scRT_pyro_supp_g_output.csv.gz'
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pert'
    log: 'logs/sig_tumors/{dataset}/infer_scRT_pyro.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/infer_scRT.py '
        '{input} {params} {output} &> {log}'


rule plot_pyro_model_output_st:
    input:
        s_phase = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT.csv.gz',
        g1_phase = 'analysis/sig_tumors/{dataset}/g1_phase_cells_with_scRT.csv.gz'
    output: 'plots/sig_tumors/{dataset}/inferred_cn_rep_results.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_pyro_model_output.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log}'


rule predict_cycle_phase_st:
    input: 
        cn_s = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT.csv.gz',
        cn_g = 'analysis/sig_tumors/{dataset}/g1_phase_cells_with_scRT.csv.gz'
    output: 
        out_s = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        out_g = 'analysis/sig_tumors/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz',
        out_lowqual = 'analysis/sig_tumors/{dataset}/model_lowqual_cells.csv.gz'
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/sig_tumors/{dataset}/predict_cycle_phase.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/predict_cycle_phase.py '
        '{input} {params} {output} &> {log}'


rule plot_filtered_pyro_model_output_st:
    input:
        s_phase = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        g1_phase = 'analysis/sig_tumors/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz'
    output: 'plots/sig_tumors/{dataset}/inferred_cn_rep_results_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_filtered_pyro_model_output.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log}'


rule compute_rt_pseudobulks_st:
    input: 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz'
    output: 'analysis/sig_tumors/{dataset}/scRT_pseudobulks.csv.gz'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/sig_tumors/{dataset}/compute_rt_pseudobulks.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log}'


rule twidth_analysis_st:
    input: 
        scrt = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        bulks = 'analysis/sig_tumors/{dataset}/scRT_pseudobulks.csv.gz'
    output: 
        output_csv = 'analysis/sig_tumors/{dataset}/twidth_values.csv.gz',
        output_png = 'plots/sig_tumors/{dataset}/twidth_curves.png'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
    log: 'logs/sig_tumors/{dataset}/twidth_analysis.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/twidth_analysis.py '
        '{input} {params} {output} &> {log}'


rule compute_cell_cycle_clone_counts_st:
    input:
        cn_s = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        cn_g = 'analysis/sig_tumors/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz'
    output: 'analysis/sig_tumors/{dataset}/cell_cycle_clone_counts.csv.gz'
    log: 'logs/sig_tumors/{dataset}/compute_cell_cycle_clone_counts.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/compute_cell_cycle_clone_counts.py '
        '{input} {output} &> {log}'


rule plot_clone_spf_st:
    input: 'analysis/sig_tumors/{dataset}/cell_cycle_clone_counts.csv.gz'
    output: 'plots/sig_tumors/{dataset}/clone_spf.png'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_clone_spf.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/plot_clone_spf.py '
        '{input} {params} {output} &> {log}'


rule plot_clone_rt_pseudobulks_st:
    input: 
        cn_s = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        rt = 'analysis/sig_tumors/{dataset}/scRT_pseudobulks.csv.gz',
        counts = 'analysis/sig_tumors/{dataset}/cell_cycle_clone_counts.csv.gz'
    output: 'plots/sig_tumors/{dataset}/clone_rt.png'
    params:
        rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_clone_rt_pseudobulks.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/plot_clone_rt_pseudobulks.py '
        '{input} {params} {output} &> {log}'


rule cohort_clone_counts_st:
    input:
        expand(
            'analysis/sig_tumors/{dataset}/cell_cycle_clone_counts.csv.gz',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
    output: 'analysis/sig_tumors/cohort_clone_counts.csv.gz'
    log: 'logs/sig_tumors/cohort_clone_counts.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/cohort_clone_counts.py '
        '--input {input} '
        '--output {output} '
        '&> {log}'


rule plot_cohort_spf_st:
    input: 'analysis/sig_tumors/cohort_clone_counts.csv.gz'
    output: 'plots/sig_tumors/cohort_spf.png'
    log: 'logs/sig_tumors/plot_cohort_spf.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/sig_tumors/plot_cohort_spf.py '
        '{input} {output} &> {log}'