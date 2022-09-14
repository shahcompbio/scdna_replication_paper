import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []
# ['D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'P1', 'P2', 'P3', 'P4', 'P5']

rule all_fig2:
    input:
        expand(
            'plots/fig2/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/scRT_heatmaps_bulk.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/twidth_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/twidth_heatmaps_bulk.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/true_scRT_heatmap.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/cn_vs_scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        


rule simulate_cell_cn_states:
    input:
        gc_rt_data = 'data/gc_rt_bin_sizes.csv',
        gc_map_data = 'data/gc_map_500kb.csv'
    output:
        s_phase = 'analysis/fig2/{dataset}/s_phase_cn_states.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cn_states.tsv'
    params:
        num_cells_S = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['num_cells_S'],
        num_cells_G = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['num_cells_G'],
        bin_size = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['bin_size'],
        clones = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['clones'],
        clone_probs = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['clone_probs'],
        states = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['states'],
        state_probs = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['state_probs'],
        cell_CNA_prob = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['cell_CNA_prob'],
        rt_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['rt_col']
    log:
        'logs/fig2/{dataset}/simulate_cell_cn_states.log'
    shell:
        'python3 scripts/fig2/simulate_cell_cn_states.py '
        '-gr {input.gc_rt_data} '
        '-gm {input.gc_map_data} '
        '-nS {params.num_cells_S} '
        '-nG {params.num_cells_G} '
        '-bs {params.bin_size} '
        '-c {params.clones} '
        '-cp {params.clone_probs} '
        '-s {params.states} '
        '-sp {params.state_probs} '
        '-cna {params.cell_CNA_prob} '
        '-rt {params.rt_col} '
        '-so {output.s_phase} '
        '-go {output.g1_phase} '
        '&> {log}'


rule simulate_reads_from_cn_pyro:
    input:
        s_phase = 'analysis/fig2/{dataset}/s_phase_cn_states.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cn_states.tsv'
    output:
        s_phase = 'analysis/fig2/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cells.tsv'
    params:
        nb_r = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['nb_r'],
        gc_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_col'],
        gc_betas = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_betas'],
        rt_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['rt_col'],
        A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
        num_reads = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['num_reads']
    log:
        'logs/fig2/{dataset}/simulate_reads_from_cn_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/simulate_reads_from_cn_pyro.py '
        '-si {input.s_phase} '
        '-gi {input.g1_phase} '
        '-nbr {params.nb_r} '
        '-gc {params.gc_col} '
        '-b {params.gc_betas} '
        '-rt {params.rt_col} '
        '-a {params.A} '
        '-n {params.num_reads} '
        '-so {output.s_phase} '
        '-go {output.g1_phase} '
        '&> {log} ; '
        'deactivate'


rule plot_cn_heatmaps:
    input:
        s_phase = 'analysis/fig2/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cells.tsv'
    output:
        s_phase = 'plots/fig2/{dataset}/cn_heatmaps.png',
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


rule plot_true_scRT_heatmap:
    input: 'analysis/fig2/{dataset}/s_phase_cells.tsv',
    output: 'plots/fig2/{dataset}/true_scRT_heatmap.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log:
        'logs/fig2/{dataset}/plot_true_scRT_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/plot_true_scRT_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule compute_ccc_features:
    input: 
        s_phase = 'analysis/fig2/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cells.tsv'
    output: 
        s_phase = 'analysis/fig2/{dataset}/s_phase_cells_features.tsv',
        g1_phase = 'analysis/fig2/{dataset}/g1_phase_cells_features.tsv'
    log: 'logs/fig2/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_bulk:
    input:
        cn_s = 'analysis/fig2/{dataset}/s_phase_cells_features.tsv',
        cn_g1 = 'analysis/fig2/{dataset}/g1_phase_cells_features.tsv'
    output: 'analysis/fig2/{dataset}/s_phase_cells_bulk_infered.tsv',
    params:
        input_col = 'true_reads_norm',
        cn_col = 'observed_cn_state',
        gc_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_col'],
        rt_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['rt_col'],
        infer_mode = 'bulk'
    log: 'logs/fig2/{dataset}/infer_scRT_bulk.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro:
    input:
        cn_s = 'analysis/fig2/{dataset}/s_phase_cells_features.tsv',
        cn_g1 = 'analysis/fig2/{dataset}/g1_phase_cells_features.tsv'
    output: 'analysis/fig2/{dataset}/s_phase_cells_pyro_infered.tsv',
    params:
        input_col = 'true_reads_norm',
        cn_col = 'observed_cn_state',
        gc_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_col'],
        rt_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['rt_col'],
        infer_mode = 'pyro'
    log: 'logs/fig2/{dataset}/infer_scRT_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule remove_nonreplicating_cells:
    input: 
        cn_pyro = 'analysis/fig2/{dataset}/s_phase_cells_pyro_infered.tsv',
        cn_bulk = 'analysis/fig2/{dataset}/s_phase_cells_bulk_infered.tsv'
    output:
        cn_pyro = 'analysis/fig2/{dataset}/s_phase_cells_pyro_filtered.tsv',
        cn_bulk = 'analysis/fig2/{dataset}/s_phase_cells_bulk_filtered.tsv'
    params:
        frac_rt_col = 'cell_frac_rep',
        pyro_rep_col = 'model_rep_state',
        bulk_rep_col = 'rt_state'
    log: 'logs/fig2/{dataset}/remove_nonreplicating_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/remove_nonreplicating_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features:
    input:
        cn_pyro = 'analysis/fig2/{dataset}/s_phase_cells_pyro_infered.tsv',
        cn_bulk = 'analysis/fig2/{dataset}/s_phase_cells_bulk_infered.tsv',
        cn_g1 = 'analysis/fig2/{dataset}/g1_phase_cells_features.tsv'
    output: 
        plot1 = 'plots/fig2/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/fig2/{dataset}/ccc_features_scatter.png'
    params:
        frac_rt_col = 'cell_frac_rep',
        pyro_rep_col = 'model_rep_state',
        bulk_rep_col = 'rt_state'
    log: 'logs/fig2/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule evaluate_model_performance_bulk:
    input: 'analysis/fig2/{dataset}/s_phase_cells_bulk_filtered.tsv'
    output: 
        plot1 = 'plots/fig2/{dataset}/scRT_heatmaps_bulk.png',
        plot2 = 'plots/fig2/{dataset}/scRT_accuracy_heatamps_bulk.png',
        plot3 = 'plots/fig2/{dataset}/frac_rt_distributions_bulk.png'
    params:
        rep_col = 'rt_state',
        cn_col = 'changepoint_segments',
        frac_rt_col = 'frac_rt'
    log: 'logs/fig2/{dataset}/evaluate_model_performance_bulk.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/evaluate_model_performance.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule evaluate_model_performance_pyro:
    input: 'analysis/fig2/{dataset}/s_phase_cells_pyro_filtered.tsv'
    output: 
        plot1 = 'plots/fig2/{dataset}/scRT_heatmaps_pyro.png',
        plot2 = 'plots/fig2/{dataset}/scRT_accuracy_heatamps_pyro.png',
        plot3 = 'plots/fig2/{dataset}/frac_rt_distributions_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/fig2/{dataset}/evaluate_model_performance_pyro.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/evaluate_model_performance.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_pyro_inferred_cn_vs_scRT:
    input: 'analysis/fig2/{dataset}/s_phase_cells_pyro_filtered.tsv'
    output: 'plots/fig2/{dataset}/cn_vs_scRT_heatmaps_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/fig2/{dataset}/plot_pyro_inferred_cn_vs_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/plot_pyro_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'



rule compute_rt_pseudobulks_pyro:
    input: 'analysis/fig2/{dataset}/s_phase_cells_pyro_filtered.tsv'
    output: 'analysis/fig2/{dataset}/scRT_pseudobulks_pyro.tsv'
    params:
        rep_col = 'model_rep_state',
        true_rep_col = 'true_rep'
    log: 'logs/fig2/{dataset}/compute_rt_pseudobulks_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_bulk:
    input: 'analysis/fig2/{dataset}/s_phase_cells_bulk_filtered.tsv'
    output: 'analysis/fig2/{dataset}/scRT_pseudobulks_bulk.tsv'
    params:
        rep_col = 'rt_state',
        true_rep_col = 'true_rep'
    log: 'logs/fig2/{dataset}/compute_rt_pseudobulks_bulk.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig2/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_pyro:
    input: 
        cn = 'analysis/fig2/{dataset}/s_phase_cells_pyro_filtered.tsv',
        pseudobulk = 'analysis/fig2/{dataset}/scRT_pseudobulks_pyro.tsv'
    output: 
        plot1 = 'plots/fig2/{dataset}/twidth_heatmaps_pyro.png',
        plot2 = 'plots/fig2/{dataset}/twidth_curves_pyro.png',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        nb_r = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['nb_r'],
        A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
        frac_rt_col = 'cell_frac_rep',
        true_frac_col = 'true_t',
        rep_state = 'model_rep_state',
        true_rep_state = 'true_rep',
        infer_mode = 'pyro'
    log: 'logs/fig2/{dataset}/twidth_analysis_pyro.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_bulk:
    input:
        cn = 'analysis/fig2/{dataset}/s_phase_cells_bulk_filtered.tsv',
        pseudobulk = 'analysis/fig2/{dataset}/scRT_pseudobulks_bulk.tsv'
    output: 
        plot1 = 'plots/fig2/{dataset}/twidth_heatmaps_bulk.png',
        plot2 = 'plots/fig2/{dataset}/twidth_curves_bulk.png',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        nb_r = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['nb_r'],
        A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
        frac_rt_col = 'frac_rt',
        true_frac_col = 'true_t',
        rep_state = 'rt_state',
        true_rep_state = 'true_rep',
        infer_mode = 'bulk'
    log: 'logs/fig2/{dataset}/twidth_analysis_bulk.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig2/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'

