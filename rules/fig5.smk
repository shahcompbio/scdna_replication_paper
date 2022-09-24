import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []

rule all_fig5:
    input:
        expand(
            'plots/fig5/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        # expand(
        #     'plots/fig5/{dataset}/rt_heatmap.png',
        #     dataset=[
        #         d for d in config['signatures_patient_tumors']
        #         if (d not in bad_datasets)
        #     ]
        # ),
        expand(
            'plots/fig5/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),


rule collect_cn_data_5:
    input: 
        hmm = 'data/signatures/signatures-hmmcopy.csv',
        annotation = 'data/signatures/signatures-annotation.csv'
    output: 'analysis/fig5/{dataset}/cn_data.tsv'
    log: 'logs/fig5/{dataset}/collect_cn_data.log'
    params:
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python scripts/fig5/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--dataset {params.dataset} '
        '--output {output} &> {log}'


rule clone_assignments_5:
    input: 
        cn = 'analysis/fig5/{dataset}/cn_data.tsv',
        clones = 'data/signatures/clone_trees/{dataset}_clones.tsv'
    output: 'analysis/fig5/{dataset}/cn_data_clones.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/fig5/{dataset}/clone_assignments.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/clone_assignments.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_ccc_features_5:
    input: 'analysis/fig5/{dataset}/cn_data_clones.tsv'
    output: 'analysis/fig5/{dataset}/cn_data_features.tsv'
    log: 'logs/fig5/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_5:
    input: 'analysis/fig5/{dataset}/cn_data_features.tsv'
    output: 
        plot1 = 'plots/fig5/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/fig5/{dataset}/ccc_features_scatter.png'
    log: 'logs/fig5/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig5/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule split_cell_cycle_5:
    input: 'analysis/fig5/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/fig5/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    log: 'logs/fig5/{dataset}/split_cell_cycle.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig5/split_cell_cycle.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_5:
    input:
        cn_s = 'analysis/fig5/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fig5/{dataset}/s_phase_cells_with_scRT.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pyro'
    log: 'logs/fig5/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_5:
    input:
        s_phase = 'analysis/fig5/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    output: 'plots/fig5/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig5/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_rt_heatmap_5:
    input: 'analysis/fig5/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'plots/fig5/{dataset}/rt_heatmap.png'
    params:
        value_col = 'model_rep_state',
        sort_col = 'model_s_time',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig5/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_pyro_model_output_5:
    input:
        s_phase = 'analysis/fig5/{dataset}/s_phase_cells_with_scRT.tsv',
        g1_phase = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fig5/{dataset}/inferred_cn_rep_results.png',
        plot2 = 'plots/fig5/{dataset}/s_vs_g_hmmcopy_states.png',
        plot3 = 'plots/fig5/{dataset}/s_vs_g_rpm.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig5/{dataset}/plot_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule remove_nonreplicating_cells_5:
    input: 'analysis/fig5/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 
        good = 'analysis/fig5/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        bad = 'analysis/fig5/{dataset}/model_nonrep_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
    log: 'logs/fig5/{dataset}/remove_nonreplicating_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/remove_nonreplicating_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_filtered_pyro_model_output_5:
    input:
        s_phase = 'analysis/fig5/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fig5/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/fig5/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/fig5/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig5/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_nonrep_pyro_model_output_5:
    input:
        s_phase = 'analysis/fig5/{dataset}/model_nonrep_cells.tsv',
        g1_phase = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fig5/{dataset}/inferred_cn_rep_results_nonrep.png',
        plot2 = 'plots/fig5/{dataset}/s_vs_g_hmmcopy_states_nonrep.png',
        plot3 = 'plots/fig5/{dataset}/s_vs_g_rpm_nonrep.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig5/{dataset}/plot_nonrep_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_5:
    input: 'analysis/fig5/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/fig5/{dataset}/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/fig5/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_5:
    input: 
        scrt = 'analysis/fig5/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        bulks = 'analysis/fig5/{dataset}/scRT_pseudobulks.tsv'
    output: 
        output_tsv = 'analysis/fig5/{dataset}/twidth_values.tsv',
        output_png = 'plots/fig5/{dataset}/twidth_curves.png'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
    log: 'logs/fig5/{dataset}/twidth_analysis.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'