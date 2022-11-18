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
        # expand(
        #     'plots/sig_tumors/{dataset}/rt_heatmap.png',
        #     dataset=[
        #         d for d in config['signatures_patient_tumors']
        #         if (d not in bad_datasets)
        #     ]
        # ),
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
            'plots/sig_tumors/{dataset}/inferred_cn_rep_results_nonrep.png',
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


rule collect_cn_data_st:
    input: 
        hmm = 'data/signatures/signatures-hmmcopy.csv',
        annotation = 'data/signatures/signatures-annotation.csv'
    output: 'analysis/sig_tumors/{dataset}/cn_data.tsv'
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
        cn = 'analysis/sig_tumors/{dataset}/cn_data.tsv',
        clones = 'data/signatures/clone_trees/{dataset}_clones.tsv'
    output: 'analysis/sig_tumors/{dataset}/cn_data_clones.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/sig_tumors/{dataset}/clone_assignments.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/clone_assignments.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_ccc_features_st:
    input: 'analysis/sig_tumors/{dataset}/cn_data_clones.tsv'
    output: 'analysis/sig_tumors/{dataset}/cn_data_features.tsv'
    log: 'logs/sig_tumors/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_st:
    input: 'analysis/sig_tumors/{dataset}/cn_data_features.tsv'
    output: 
        plot1 = 'plots/sig_tumors/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/sig_tumors/{dataset}/ccc_features_scatter.png'
    log: 'logs/sig_tumors/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_tumors/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule split_cell_cycle_st:
    input: 'analysis/sig_tumors/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/sig_tumors/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/sig_tumors/{dataset}/g1_phase_cells.tsv'
    log: 'logs/sig_tumors/{dataset}/split_cell_cycle.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_tumors/split_cell_cycle.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_st:
    input:
        cn_s = 'analysis/sig_tumors/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/sig_tumors/{dataset}/g1_phase_cells.tsv'
    output:
        main_s_out = 'analysis/laks_flow/{dataset}/s_phase_cells_with_scRT.tsv',
        supp_s_out = 'analysis/laks_flow/{dataset}/scRT_pyro_supp_s_output.tsv',
        main_g_out = 'analysis/laks_flow/{dataset}/g1_phase_cells_with_scRT.tsv',
        supp_g_out = 'analysis/laks_flow/{dataset}/scRT_pyro_supp_g_output.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pyro'
    log: 'logs/sig_tumors/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_st:
    input:
        s_phase = 'analysis/sig_tumors/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/sig_tumors/{dataset}/g1_phase_cells.tsv'
    output: 'plots/sig_tumors/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/sig_lines/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_rt_heatmap_st:
    input: 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'plots/sig_tumors/{dataset}/rt_heatmap.png'
    params:
        value_col = 'model_rep_state',
        sort_col = 'model_s_time',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/sig_lines/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_pyro_model_output_st:
    input:
        s_phase = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT.tsv',
        g1_phase = 'analysis/sig_tumors/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/sig_tumors/{dataset}/inferred_cn_rep_results.png',
        plot2 = 'plots/sig_tumors/{dataset}/s_vs_g_hmmcopy_states.png',
        plot3 = 'plots/sig_tumors/{dataset}/s_vs_g_rpm.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule remove_nonreplicating_cells_st:
    input: 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 
        good = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        bad = 'analysis/sig_tumors/{dataset}/model_nonrep_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
    log: 'logs/sig_tumors/{dataset}/remove_nonreplicating_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/remove_nonreplicating_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_filtered_pyro_model_output_st:
    input:
        s_phase = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/sig_tumors/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/sig_tumors/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/sig_tumors/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/sig_tumors/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_nonrep_pyro_model_output_st:
    input:
        s_phase = 'analysis/sig_tumors/{dataset}/model_nonrep_cells.tsv',
        g1_phase = 'analysis/sig_tumors/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/sig_tumors/{dataset}/inferred_cn_rep_results_nonrep.png',
        plot2 = 'plots/sig_tumors/{dataset}/s_vs_g_hmmcopy_states_nonrep.png',
        plot3 = 'plots/sig_tumors/{dataset}/s_vs_g_rpm_nonrep.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_tumors/{dataset}/plot_nonrep_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_st:
    input: 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/sig_tumors/{dataset}/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/sig_tumors/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_st:
    input: 
        scrt = 'analysis/sig_tumors/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        bulks = 'analysis/sig_tumors/{dataset}/scRT_pseudobulks.tsv'
    output: 
        output_tsv = 'analysis/sig_tumors/{dataset}/twidth_values.tsv',
        output_png = 'plots/sig_tumors/{dataset}/twidth_curves.png'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        infer_mode = 'pyro_composite',
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
    log: 'logs/sig_tumors/{dataset}/twidth_analysis.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/sig_lines/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'