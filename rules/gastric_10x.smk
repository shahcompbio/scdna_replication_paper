import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []

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

