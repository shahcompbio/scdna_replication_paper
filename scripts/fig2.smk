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
        s_time_dist = lambda wildcards: config['simulated_datasets']['diploid'][wildcards.dataset]['s_time_dist']
    log:
        'logs/fig2/{dataset}/simulate_diploid_data.log'
    shell:
        'python3 scripts/fig2/simulate_diploid_data.py '
        '{input} {params} {output} &> {log}'