import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []

rule all_fitness:
    input:
        expand(
            'plots/fitness_lines/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/clone_spf.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'analysis/fitness_lines/{dataset}/cn_pseudobulks.tsv',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness_lines/{dataset}/clones_vs_time.png',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),

# use the model output files from sig_lines.smk to assign S-phase cells to timepoints
rule assign_timepoints_fl:
    input: 
        cn_s = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        cn_g = 'analysis/sig_lines/{dataset}/g1_phase_cells.tsv',
        times = 'data/fitness/dlp_summaries_rebuttal.csv'
    output: 
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        cn_g = 'analysis/fitness_lines/{dataset}/g1_phase_cells.tsv',
    log: 'logs/fitness_lines/{dataset}/assign_timepoints.log'
    shell:
        'python3 scripts/fitness_lines/assign_timepoints.py '
        '{input} {output} &> {log}'


rule remove_nonreplicating_cells_fl:
    input: 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 
        good = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        nonrep = 'analysis/fitness_lines/{dataset}/model_nonrep_cells.tsv',
        lowqual = 'analysis/fitness_lines/{dataset}/model_lowqual_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/fitness_lines/{dataset}/remove_nonreplicating_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/common/remove_nonreplicating_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_filtered_pyro_model_output_sl:
    input:
        s_phase = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/fitness_lines/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fitness_lines/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/fitness_lines/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/fitness_lines/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/common/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_fl:
    input: 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/fitness_lines/{dataset}/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/fitness_lines/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness_lines/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_cn_pseudobulks_fl:
    input: 'analysis/fitness_lines/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fitness_lines/{dataset}/cn_pseudobulks.tsv'
    params:
        cn_state_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/compute_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness_lines/compute_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_fl:
    input:
        s_phase = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        g1_phase = 'analysis/fitness_lines/{dataset}/g1_phase_cells.tsv'
    output: 'plots/fitness_lines/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fitness_lines/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_clone_rt_and_spf_sl:
    input: 
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        cn_g ='analysis/fitness_lines/{dataset}/g1_phase_cells.tsv',
        cn_g_recovered ='analysis/fitness_lines/{dataset}/model_nonrep_cells.tsv',
        rt = 'analysis/fitness_lines/{dataset}/scRT_pseudobulks.tsv'
    output:
        tsv = 'analysis/fitness_lines/{dataset}/cell_cycle_clone_counts.tsv',
        clone_rt = 'plots/fitness_lines/{dataset}/clone_rt.png',
        clone_spf = 'plots/fitness_lines/{dataset}/clone_spf.png'
    params:
        rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness_lines/{dataset}/plot_clone_rt_and_spf.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness_lines/plot_clone_rt_and_spf.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_clones_vs_time:
    input: 'analysis/fitness_lines/{dataset}/cell_cycle_clone_counts.tsv'
    output: 
        plot1 = 'plots/fitness_lines/{dataset}/clones_vs_time.png',
        plot2 = 'plots/fitness_lines/{dataset}/total_cells_vs_time.png'
    log: 'logs/fitness_lines/{dataset}/plot_clones_vs_time.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness_lines/plot_clones_vs_time.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'
