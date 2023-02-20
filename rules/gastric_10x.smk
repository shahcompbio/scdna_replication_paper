import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = ['KATOIII', 'SNU-668']

rule all_gastric_10x:
    input:
        expand(
            'plots/gastric_10x/{dataset}/cn_heatmaps_input.png',
            dataset=[
                d for d in config['10x_gastric_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'analysis/gastric_10x/{dataset}/s_phase_cells_with_scRT.csv.gz',
            dataset=[
                d for d in config['10x_gastric_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/gastric_10x/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['10x_gastric_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        # expand(
        #     'plots/gastric_10x/{dataset}/inferred_cn_rep_results_filtered.png',
        #     dataset=[
        #         d for d in config['10x_gastric_cell_lines']
        #         if (d not in bad_datasets)
        #     ]
        # ),


rule collect_cn_data_g10x:
    input:
        cn_data = 'data/gastric_10x/{dataset}/cnv_data.h5',
        clones = 'data/gastric_10x/{dataset}/{dataset}_cnv_meta.csv',
    output: 'analysis/gastric_10x/{dataset}/cn_data.csv.gz'
    log: 'logs/gastric_10x/{dataset}/collect_cn_data.log'
    shell: 
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/gastric_10x/collect_cn_data.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_input_g10x:
    input: 'analysis/gastric_10x/{dataset}/cn_data.csv.gz'
    output: 'plots/gastric_10x/{dataset}/cn_heatmaps_input.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/gastric_10x/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/gastric_10x/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_g10x:
    input: 'analysis/gastric_10x/{dataset}/cn_data.csv.gz'
    output:
        main_s_out = 'analysis/gastric_10x/{dataset}/s_phase_cells_with_scRT.csv.gz',
        supp_s_out = 'analysis/gastric_10x/{dataset}/scRT_pyro_supp_s_output.csv.gz',
        main_g_out = 'analysis/gastric_10x/{dataset}/g1_phase_cells_with_scRT.csv.gz',
        supp_g_out = 'analysis/gastric_10x/{dataset}/scRT_pyro_supp_g_output.csv.gz',
    params:
        input_col = 'reads',
        clone_col = 'clone_id',
        cn_col = 'state',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pyro'
    log: 'logs/gastric_10x/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/gastric_10x/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_pyro_model_output_g10x:
    input:
        s_phase = 'analysis/gastric_10x/{dataset}/s_phase_cells_with_scRT.csv.gz',
        g1_phase = 'analysis/gastric_10x/{dataset}/g1_phase_cells_with_scRT.csv.gz'
    output:
        plot1 = 'plots/gastric_10x/{dataset}/inferred_cn_rep_results.png',
        plot2 = 'plots/gastric_10x/{dataset}/s_vs_g_hmmcopy_states.png',
        plot3 = 'plots/gastric_10x/{dataset}/s_vs_g_rpm.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/gastric_10x/{dataset}/plot_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/gastric_10x/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule revise_cell_cycle_labels_g10x:
    input: 
        cn_s = 'analysis/gastric_10x/{dataset}/s_phase_cells_with_scRT.csv.gz',
        cn_g = 'analysis/gastric_10x/{dataset}/g1_phase_cells_with_scRT.csv.gz',
    output:
        out_s = 'analysis/gastric_10x/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        out_g = 'analysis/gastric_10x/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz',
        out_lowqual = 'analysis/gastric_10x/{dataset}/model_lowqual_cells.csv.gz',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/gastric_10x/{dataset}/revise_cell_cycle_labels.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/gastric_10x/revise_cell_cycle_labels.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_filtered_pyro_model_output_g10x:
    input:
        s_phase = 'analysis/gastric_10x/{dataset}/s_phase_cells_with_scRT_filtered.csv.gz',
        g1_phase = 'analysis/gastric_10x/{dataset}/g1_phase_cells_with_scRT_filtered.csv.gz'
    output:
        plot1 = 'plots/gastric_10x/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/gastric_10x/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/gastric_10x/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/gastric_10x/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/gastric_10x/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'
        