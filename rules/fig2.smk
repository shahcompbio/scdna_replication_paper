import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []

rule all_fig2:
    input:
        expand(
            'analysis/fig2/{dataset}/s_phase_cells.tsv',
            dataset=[
                d for d in config['simulated_datasets']['diploid']
                if (d not in bad_datasets)
            ]
        ),
        # expand(
        #     'plots/fig2/{dataset}/scRT_heatmaps.pdf',
        #     dataset=[
        #         d for d in config['simulated_datasets']['diploid']
        #         if (d not in bad_datasets)
        #     ]
        # ),
        # expand(
        #     'plots/fig2/{dataset}/twidth_heatmaps.pdf',
        #     dataset=[
        #         d for d in config['simulated_datasets']['diploid']
        #         if (d not in bad_datasets)
        #     ]
        # ),
        expand(
            'plots/fig2/{dataset}/cn_heatmaps.pdf',
            dataset=[
                d for d in config['simulated_datasets']['diploid']
                if (d not in bad_datasets)
            ]
        ),

rule simulate_diploid_data:
    input:
        gc_rt_data = 'data/gc_rt_bin_sizes.csv',
        gc_map_data = 'data/gc_map_500kb.csv'
    output:
        s_phase = 'analysis/fig2/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cells.tsv'
    params:
        sigma1 = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['sigma1'],
        gc_slope = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['gc_slope'],
        gc_int = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['gc_int'],
        A = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['A'],
        B = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['B'],
        num_reads = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['num_reads'],
        num_cells_S = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['num_cells_S'],
        num_cells_G = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['num_cells_G'],
        bin_size = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['bin_size'],
        s_time_stdev = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['s_time_stdev']
    log:
        'logs/fig2/{dataset}/simulate_diploid_data.log'
    shell:
        'python3 scripts/fig2/simulate_diploid_data.py '
        '{input} {params} {output} &> {log}'


rule plot_cn_heatmaps:
    input:
        s_phase = 'analysis/fig2/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cells.tsv'
    output:
        s_phase = 'plots/fig2/{dataset}/cn_heatmaps.pdf',
    params:
        value_col = 'true_G1_state',
        dataset = lambda wildcards: wildcards.dataset
    log:
        'logs/fig2/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule infer_scRT:
    input:
        cn_s = 'analysis/fig2/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig2/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fig2/{dataset}/s_phase_cells_with_scRT.tsv',
    params:
        input_col = 'reads'
    log: 'logs/fig2/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule evaluate_model_performance:
    input: 'analysis/fig2/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 
        plot1 = 'plots/fig2/{dataset}/scRT_heatmaps.pdf',
        plot2 = 'plots/fig2/{dataset}/scRT_accuracy_heatamps.pdf',
        plot3 = 'plots/fig2/{dataset}/frac_rt_distributions.pdf'
    log: 'logs/fig2/{dataset}/evaluate_model_performance.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/evaluate_model_performance.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis:
    input: 'analysis/fig2/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 
        plot1 = 'plots/fig2/{dataset}/twidth_heatmaps.pdf',
        plot2 = 'plots/fig2/{dataset}/twidth_curves.pdf',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        sigma1 = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['sigma1'],
        gc_slope = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['gc_slope'],
        gc_int = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['gc_int'],
        A = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['A'],
        s_time_stdev = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['s_time_stdev']
    log: 'logs/fig2/{dataset}/twidth_analysis.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'

