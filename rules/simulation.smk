import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []
datasets = ['D1.0']
# ['D2', 'D3', 'D4', 'D5', 'D6', 'D7', 'D8', 'P1', 'P2', 'P3', 'P4', 'P5']

rule all_simulation:
    input:
        expand(
            'plots/simulation/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/scRT_heatmaps_kronos.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/scRT_heatmaps_pyro_composite.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/twidth_curves_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/twidth_curves_kronos.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/twidth_curves_pyro_composite.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/true_scRT_heatmap.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/cn_vs_scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/cn_vs_scRT_composite_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        'plots/simulation/P5.8/true_vs_inferred_heatmaps.png',
        'plots/simulation/all/model_accuracies1.png',
        'plots/simulation/all/clone_specific_rt_corr.png',
        

rule simulate_cell_cn_states_sim:
    input:
        gc_rt_data = 'data/gc_rt_bins.csv',
        gc_map_data = 'data/gc_map_500kb.csv'
    output:
        s_phase = 'analysis/simulation/{dataset}/s_phase_cn_states.tsv',
        g1_phase = 'analysis/simulation/{dataset}/g1_phase_cn_states.tsv'
    params:
        num_cells_S = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['num_cells_S'],
        num_cells_G = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['num_cells_G'],
        bin_size = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['bin_size'],
        clones = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['clones'],
        clone_probs = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['clone_probs'],
        states = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['states'],
        state_probs = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['state_probs'],
        cell_CNA_prob = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['cell_CNA_prob'],
    log:
        'logs/simulation/{dataset}/simulate_cell_cn_states.log'
    shell:
        'python3 scripts/simulation/simulate_cell_cn_states.py '
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
        '-so {output.s_phase} '
        '-go {output.g1_phase} '
        '&> {log}'


rule simulate_reads_from_cn_pyro_sim:
    input:
        s_phase = 'analysis/simulation/{dataset}/s_phase_cn_states.tsv',
        g1_phase = 'analysis/simulation/{dataset}/g1_phase_cn_states.tsv'
    output:
        s_phase = 'analysis/simulation/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/simulation/{dataset}/g1_phase_cells.tsv'
    params:
        lamb = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['lambda'],
        gc_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_col'],
        gc_betas = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_betas'],
        rt_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['rt_col'],
        A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
        num_reads = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['num_reads'],
        clones = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['clones']
    log:
        'logs/simulation/{dataset}/simulate_reads_from_cn_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/simulate_reads_from_cn_pyro.py '
        '-si {input.s_phase} '
        '-gi {input.g1_phase} '
        '-l {params.lamb} '
        '-gc {params.gc_col} '
        '-b {params.gc_betas} '
        '-rt {params.rt_col} '
        '-a {params.A} '
        '-n {params.num_reads} '
        '-c {params.clones} '
        '-so {output.s_phase} '
        '-go {output.g1_phase} '
        '&> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_sim:
    input:
        s_phase = 'analysis/simulation/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/simulation/{dataset}/g1_phase_cells.tsv'
    output:
        s_phase = 'plots/simulation/{dataset}/cn_heatmaps.png',
    params:
        value_col = 'true_G1_state',
        dataset = lambda wildcards: wildcards.dataset
    log:
        'logs/simulation/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/simulation/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_true_scRT_heatmap_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells.tsv',
    output: 'plots/simulation/{dataset}/true_scRT_heatmap.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log:
        'logs/simulation/{dataset}/plot_true_scRT_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/simulation/plot_true_scRT_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule compute_ccc_features_sim:
    input: 
        s_phase = 'analysis/simulation/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/simulation/{dataset}/g1_phase_cells.tsv'
    output: 
        s_phase = 'analysis/simulation/{dataset}/s_phase_cells_features.tsv',
        g1_phase = 'analysis/simulation/{dataset}/g1_phase_cells_features.tsv'
    log: 'logs/simulation/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# rule infer_scRT_bulk_sim:
#     input:
#         cn_s = 'analysis/simulation/{dataset}/s_phase_cells_features.tsv',
#         cn_g1 = 'analysis/simulation/{dataset}/g1_phase_cells_features.tsv'
#     output: 
#         main_s_out = 'analysis/simulation/{dataset}/s_phase_cells_bulk_inferred.tsv',
#         supp_s_out = 'analysis/simulation/{dataset}/scRT_bulk_supp_s_output.tsv',
#         main_g_out = 'analysis/simulation/{dataset}/g1_phase_cells_bulk_inferred.tsv',
#         supp_g_out = 'analysis/simulation/{dataset}/scRT_bulk_supp_g_output.tsv',
#     params:
#         input_col = 'true_reads_norm',
#         cn_col = 'observed_cn_state',
#         gc_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_col'],
#         cn_prior_method = 'diploid',  # irrelevant bc pyro model not invoked here
#         infer_mode = 'bulk',
#         max_iter = 2000
#     log: 'logs/simulation/{dataset}/infer_scRT_bulk.log'
#     shell:
#         'source ../scdna_replication_tools/venv3/bin/activate ; '
#         'python3 scripts/simulation/infer_scRT.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'


rule infer_kronos_scRT_sim:
    input:
        cn_s = 'analysis/simulation/{dataset}/s_phase_cells_features.tsv',
        cn_g1 = 'analysis/simulation/{dataset}/g1_phase_cells_features.tsv'
    output: 'analysis/simulation/{dataset}/s_phase_cells_kronos_output.tsv'
    conda: '../envs/Kronos_scRT.yaml'
    log: 'logs/simulation/{dataset}/infer_kronos_scRT.log'
    shell:
        'Kronos RT '
        '-s {input.cn_s} '
        '-g {input.cn_g1} '
        '-o {output} '
        '&> {log}'


rule process_kronos_output_sim:
    input: 
        kronos_input = 'analysis/simulation/{dataset}/s_phase_cells_features.tsv',
        kronos_output = 'analysis/simulation/{dataset}/s_phase_cells_kronos_output.tsv'
    output: 'analysis/simulation/{dataset}/s_phase_cells_kronos_inferred.tsv'
    log: 'logs/simulation/{dataset}/process_kronos_output.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/process_kronos_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_sim:
    input:
        cn_s = 'analysis/simulation/{dataset}/s_phase_cells_features.tsv',
        cn_g1 = 'analysis/simulation/{dataset}/g1_phase_cells_features.tsv'
    output:
        main_s_out = 'analysis/simulation/{dataset}/s_phase_cells_pyro_inferred.tsv',
        supp_s_out = 'analysis/simulation/{dataset}/scRT_pyro_supp_s_output.tsv',
        main_g_out = 'analysis/simulation/{dataset}/g1_phase_cells_pyro_inferred.tsv',
        supp_g_out = 'analysis/simulation/{dataset}/scRT_pyro_supp_g_output.tsv',
    params:
        input_col = 'true_reads_norm',
        cn_col = 'observed_cn_state',
        gc_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_col'],
        cn_prior_method = 'g1_clones',
        infer_mode = 'pyro',
        max_iter = 2000
    log: 'logs/simulation/{dataset}/infer_scRT_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_composite_sim:
    input:
        cn_s = 'analysis/simulation/{dataset}/s_phase_cells_features.tsv',
        cn_g1 = 'analysis/simulation/{dataset}/g1_phase_cells_features.tsv'
    output:
        main_s_out = 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv',
        supp_s_out = 'analysis/simulation/{dataset}/scRT_pyro_composite_supp_s_output.tsv',
        main_g_out = 'analysis/simulation/{dataset}/g1_phase_cells_pyro_composite_inferred.tsv',
        supp_g_out = 'analysis/simulation/{dataset}/scRT_pyro_composite_supp_g_output.tsv',
    params:
        input_col = 'true_reads_norm',
        cn_col = 'observed_cn_state',
        gc_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['gc_col'],
        cn_prior_method = 'g1_composite',
        infer_mode = 'pyro',
        max_iter = 2000
    log: 'logs/simulation/{dataset}/infer_scRT_pyro_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# rule remove_nonreplicating_cells_sim:
#     input: 
#         cn_pyro = 'analysis/simulation/{dataset}/s_phase_cells_pyro_inferred.tsv',
#         cn_kronos = 'analysis/simulation/{dataset}/s_phase_cells_kronos_inferred.tsv'
#     output:
#         cn_pyro = 'analysis/simulation/{dataset}/s_phase_cells_pyro_filtered.tsv',
#         cn_kronos = 'analysis/simulation/{dataset}/s_phase_cells_kronos_filtered.tsv'
#     params:
#         frac_rt_col = 'cell_frac_rep',
#         pyro_rep_col = 'model_rep_state',
#         kronos_rep_col = 'rt_state'
#     log: 'logs/simulation/{dataset}/remove_nonreplicating_cells.log'
#     shell:
#         'source ../scdna_replication_tools/venv3/bin/activate ; '
#         'python3 scripts/simulation/remove_nonreplicating_cells.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'


# rule remove_nonreplicating_cells_composite_sim:
#     input: 
#         cn_pyro = 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv',
#         cn_kronos = 'analysis/simulation/{dataset}/s_phase_cells_kronos_inferred.tsv'
#     output:
#         cn_pyro = 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_filtered.tsv',
#         cn_kronos = 'analysis/simulation/{dataset}/s_phase_cells_kronos_composite_filtered.tsv'
#     params:
#         frac_rt_col = 'cell_frac_rep',
#         pyro_rep_col = 'model_rep_state',
#         kronos_rep_col = 'rt_state'
#     log: 'logs/simulation/{dataset}/remove_nonreplicating_cells_composite.log'
#     shell:
#         'source ../scdna_replication_tools/venv3/bin/activate ; '
#         'python3 scripts/simulation/remove_nonreplicating_cells.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'


rule revise_cell_cycle_labels_composite_sim:
    input: 
        cn_s = 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv',
        cn_g = 'analysis/simulation/{dataset}/g1_phase_cells_pyro_composite_inferred.tsv',
    output:
        out_s = 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_filtered.tsv',
        out_g = 'analysis/simulation/{dataset}/g1_phase_cells_pyro_composite_filtered.tsv',
        out_lowqual = 'analysis/simulation/{dataset}/model_lowqual_composite_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'true_reads_norm'
    log: 'logs/simulation/{dataset}/revise_cell_cycle_labels.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/revise_cell_cycle_labels.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_sim:
    input:
        cn_s = 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_filtered.tsv',
        cn_g1 = 'analysis/simulation/{dataset}/g1_phase_cells_pyro_composite_filtered.tsv',
        cn_lowqual = 'analysis/simulation/{dataset}/model_lowqual_composite_cells.tsv'
    output: 
        plot1 = 'plots/simulation/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/simulation/{dataset}/ccc_features_scatter.png',
        plot3 = 'plots/simulation/{dataset}/predicted_phase_confusion_mat.png'
    params:
        frac_rt_col = 'cell_frac_rep',
        pyro_rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/simulation/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule evaluate_model_performance_kronos_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_kronos_inferred.tsv'
    output: 
        plot1 = 'plots/simulation/{dataset}/scRT_heatmaps_kronos.png',
        plot2 = 'plots/simulation/{dataset}/scrt_accuracy_heatmaps_kronos.png',
        plot3 = 'plots/simulation/{dataset}/frac_rt_distributions_kronos.png'
    params:
        rep_col = 'rt_state',
        cn_col = 'observed_cn_state',
        frac_rt_col = 'frac_rt'
    log: 'logs/simulation/{dataset}/evaluate_model_performance_kronos.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/simulation/evaluate_model_performance.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule evaluate_model_performance_pyro_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_pyro_inferred.tsv'
    output: 
        plot1 = 'plots/simulation/{dataset}/scRT_heatmaps_pyro.png',
        plot2 = 'plots/simulation/{dataset}/scrt_accuracy_heatmaps_pyro.png',
        plot3 = 'plots/simulation/{dataset}/frac_rt_distributions_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/simulation/{dataset}/evaluate_model_performance_pyro.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/simulation/evaluate_model_performance.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule evaluate_model_performance_pyro_composite_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv'
    output: 
        plot1 = 'plots/simulation/{dataset}/scRT_heatmaps_pyro_composite.png',
        plot2 = 'plots/simulation/{dataset}/scrt_accuracy_heatmaps_pyro_composite.png',
        plot3 = 'plots/simulation/{dataset}/frac_rt_distributions_pyro_composite.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/simulation/{dataset}/evaluate_model_performance_pyro_composite.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/simulation/evaluate_model_performance.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_pyro_inferred_cn_vs_scRT_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_pyro_inferred.tsv'
    output: 'plots/simulation/{dataset}/cn_vs_scRT_heatmaps_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/simulation/{dataset}/plot_pyro_inferred_cn_vs_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/plot_pyro_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_pyro_composite_inferred_cn_vs_scRT_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv'
    output: 'plots/simulation/{dataset}/cn_vs_scRT_composite_heatmaps_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/simulation/{dataset}/plot_pyro_composite_inferred_cn_vs_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/plot_pyro_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule true_vs_inferred_heatmaps_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv',
    output: 'plots/simulation/{dataset}/true_vs_inferred_heatmaps.png',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        true_frac_col = 'true_t',
        rep_state = 'model_rep_state',
        true_rep_state = 'true_rep'
    log: 'logs/simulation/{dataset}/true_vs_inferred_heatmaps.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/true_vs_inferred_heatmaps.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_pyro_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_pyro_inferred.tsv'
    output: 'analysis/simulation/{dataset}/scRT_pseudobulks_pyro.tsv'
    params:
        rep_col = 'model_rep_state',
        true_rep_col = 'true_rep'
    log: 'logs/simulation/{dataset}/compute_rt_pseudobulks_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_kronos_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_kronos_inferred.tsv'
    output: 'analysis/simulation/{dataset}/scRT_pseudobulks_kronos.tsv'
    params:
        rep_col = 'rt_state',
        true_rep_col = 'true_rep'
    log: 'logs/simulation/{dataset}/compute_rt_pseudobulks_kronos.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_pyro_composite_sim:
    input: 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv'
    output: 'analysis/simulation/{dataset}/scRT_pseudobulks_pyro_composite.tsv'
    params:
        rep_col = 'model_rep_state',
        true_rep_col = 'true_rep'
    log: 'logs/simulation/{dataset}/compute_rt_pseudobulks_pyro_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_pyro_sim:
    input: 
        cn = 'analysis/simulation/{dataset}/s_phase_cells_pyro_inferred.tsv',
        pseudobulk = 'analysis/simulation/{dataset}/scRT_pseudobulks_pyro.tsv'
    output: 
        tsv = 'analysis/simulation/{dataset}/twidth_values_pyro.tsv',
        plot = 'plots/simulation/{dataset}/twidth_curves_pyro.png',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        lamb = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['lambda'],
        A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
        frac_rt_col = 'cell_frac_rep',
        true_frac_col = 'true_t',
        rep_state = 'model_rep_state',
        true_rep_state = 'true_rep',
        infer_mode = 'pyro'
    log: 'logs/simulation/{dataset}/twidth_analysis_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_kronos_sim:
    input:
        cn = 'analysis/simulation/{dataset}/s_phase_cells_kronos_inferred.tsv',
        pseudobulk = 'analysis/simulation/{dataset}/scRT_pseudobulks_kronos.tsv'
    output: 
        tsv = 'analysis/simulation/{dataset}/twidth_values_kronos.tsv',
        plot = 'plots/simulation/{dataset}/twidth_curves_kronos.png',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        lamb = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['lambda'],
        A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
        frac_rt_col = 'frac_rt',
        true_frac_col = 'true_t',
        rep_state = 'rt_state',
        true_rep_state = 'true_rep',
        infer_mode = 'kronos'
    log: 'logs/simulation/{dataset}/twidth_analysis_kronos.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_pyro_composite_sim:
    input: 
        cn = 'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_inferred.tsv',
        pseudobulk = 'analysis/simulation/{dataset}/scRT_pseudobulks_pyro_composite.tsv'
    output: 
        tsv = 'analysis/simulation/{dataset}/twidth_values_pyro_composite.tsv',
        plot = 'plots/simulation/{dataset}/twidth_curves_pyro_composite.png',
    params:
        dataset = lambda wildcards: wildcards.dataset,
        lamb = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['lambda'],
        A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
        frac_rt_col = 'cell_frac_rep',
        true_frac_col = 'true_t',
        rep_state = 'model_rep_state',
        true_rep_state = 'true_rep',
        infer_mode = 'pyro_composite'
    log: 'logs/simulation/{dataset}/twidth_analysis_pyro_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# takeaway is that rule_aggregate_model_results succeeds when =<20 files are provided 
# but fails when >20 files are provided
# dataset_list = [d for d in config['simulated_datasets'] if d.startswith('D1.0')]
# dataset_list = ['D1.0', 'D1.1', 'D1.2', 'D1.3', 'D1.4', 'D1.5', 'D1.6', 'D1.7', 'D1.8']
# dataset_list = ['D2.0', 'D2.1', 'D2.2', 'D2.3', 'D2.4', 'D2.5']

# def get_cn_bulk_paths():
#     files = expand(
#         'analysis/simulation/{dataset}/s_phase_cells_bulk_filtered.tsv',
#         dataset=dataset_list
#     )
#     return files


# def get_cn_pyro_clone_paths():
#     files = expand(
#         'analysis/simulation/{dataset}/s_phase_cells_pyro_filtered.tsv',
#         dataset=dataset_list
#     )
#     return files


# def get_cn_pyro_comp_paths():
#     files = expand(
#         'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_filtered.tsv',
#         dataset=dataset_list
#     )
#     return files


# dataset_list = ['D1.6', 'D1.7', 'D1.8']
# print(dataset_list)
# rule aggregate_model_results:
#     input: 
#         cn_bulk = expand(
#             'analysis/simulation/{dataset}/s_phase_cells_bulk_filtered.tsv',
#             dataset=dataset_list
#         ),
#         cn_pyro_clone = expand(
#             'analysis/simulation/{dataset}/s_phase_cells_pyro_filtered.tsv',
#             dataset=dataset_list
#         ),
#         cn_pyro_comp = expand(
#             'analysis/simulation/{dataset}/s_phase_cells_pyro_composite_filtered.tsv',
#             dataset=dataset_list
#         )
#     output: 'analysis/simulation/all/s_phase_model_results_paths.tsv'
#     params:
#         datasets = expand(dataset_list)
#     run:
#         df = []
#         for dataset, bulk_path, clone_path, comp_path in zip(
#                 params.datasets, input.cn_bulk, input.cn_pyro_clone, input.cn_pyro_comp
#             ):
#             temp_df = pd.DataFrame({
#                 'dataset': [str(dataset)], 'bulk_path': [str(bulk_path)],
#                 'clone_path': [str(clone_path)], 'comp_path': [str(comp_path)]
#             })
#             df.append(temp_df)
#         df = pd.concat(df, ignore_index=True)
#         df.to_csv(str(output), sep='\t', index=False)


rule aggregate_model_results_sim:
    input: 'analysis/simulation/D1.0/s_phase_cells_pyro_composite_inferred.tsv'
    output: 'analysis/simulation/all/s_phase_model_results_paths.tsv'
    params:
        datasets = expand(config['simulated_datasets'])
    log: 'logs/simulation/all/aggregate_model_results.log'
    shell:
        'python3 scripts/simulation/aggregate_model_results.py '
        '-d {params.datasets} '
        '-o {output} '
         '&> {log}'


rule model_accuracies_sim:
    input: 'analysis/simulation/all/s_phase_model_results_paths.tsv'
    output:
        accuracy_table = 'analysis/simulation/all/model_accuracies.tsv',
        accuracy_plot1 = 'plots/simulation/all/model_accuracies1.png',
        accuracy_plot2 = 'plots/simulation/all/model_accuracies2.png',
    params:
        datasets = expand(config['simulated_datasets']),
        A = expand([str(config['simulated_datasets'][d]['A']) for d in config['simulated_datasets']]),
        cell_cna_rate = expand([str(config['simulated_datasets'][d]['cell_CNA_prob']) for d in config['simulated_datasets']]),
        num_clones = expand([str(len(config['simulated_datasets'][d]['clones'])) for d in config['simulated_datasets']]),
        lamb = expand([str(config['simulated_datasets'][d]['lambda']) for d in config['simulated_datasets']]),
        kronos_rep_col = 'rt_state',
        pyro_rep_col = 'model_rep_state',
        pyro_cn_col = 'model_cn_state',
        true_rep_col = 'true_rep',
        true_cn_col = 'true_G1_state'
    log: 'logs/simulation/all/model_accuracies.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/model_accuracies.py '
        '--input {input} '
        '--datasets {params.datasets} '
        '--A {params.A} '
        '--cell_cna_rate {params.cell_cna_rate} '
        '--num_clones {params.num_clones} '
        '--lamb {params.lamb} '
        '--kronos_rep_col {params.kronos_rep_col} '
        '--pyro_rep_col {params.pyro_rep_col} '
        '--pyro_cn_col {params.pyro_cn_col} '
        '--true_rep_col {params.true_rep_col} '
        '--true_cn_col {params.true_cn_col} '
        '--table {output.accuracy_table} '
        '--plot1 {output.accuracy_plot1} '
        '--plot2 {output.accuracy_plot2} '
        '&> {log} ; '
        'deactivate'


rule clone_specific_rt_sim:
    input:
        P10 = 'analysis/simulation/P10/scRT_pseudobulks_pyro_composite.tsv',
        P11 = 'analysis/simulation/P11/scRT_pseudobulks_pyro_composite.tsv'
    output: 
        corr = 'plots/simulation/all/clone_specific_rt_corr.png',
        profiles = 'plots/simulation/all/clone_specific_rt_profiles.png'
    log: 'logs/simulation/all/clone_specific_rt.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/simulation/clone_specific_rt.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'
