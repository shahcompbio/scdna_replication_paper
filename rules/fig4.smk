import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

datasets = ['all', 'GM18507', 'T47D']

training_data = pd.read_csv('/juno/work/shah/users/weinera2/projects/cell_cycle_classifier/cell_cycle_classifier/data/training/curated_feature_data_rt.csv')

# paths on juno for reads, metics, and annotation metrics
cn_data_urls = [
    '/work/shah/tantalus/SC-1563/results/results/hmmcopy_autoploidy/A90553C_multiplier0_reads.csv.gz',
    '/work/shah/tantalus/SC-1561/results/results/hmmcopy_autoploidy/A73044A_multiplier0_reads.csv.gz',
    '/work/shah/tantalus/SC-1583/results/results/hmmcopy_autoploidy/A96139A_multiplier0_reads.csv.gz',
    '/work/shah/tantalus/SC-1585/results/results/hmmcopy_autoploidy/A96147A_multiplier0_reads.csv.gz',
]

metrics_data_urls = [
    '/work/shah/tantalus/SC-1563/results/results/hmmcopy_autoploidy/A90553C_multiplier0_metrics.csv.gz',
    '/work/shah/tantalus/SC-1561/results/results/hmmcopy_autoploidy/A73044A_multiplier0_metrics.csv.gz',
    '/work/shah/tantalus/SC-1583/results/results/hmmcopy_autoploidy/A96139A_multiplier0_metrics.csv.gz',
    '/work/shah/tantalus/SC-1585/results/results/hmmcopy_autoploidy/A96147A_multiplier0_metrics.csv.gz',
]

align_metrics_data_urls = [
    '/work/shah/tantalus/SC-1563/results/results/alignment/A90553C_alignment_metrics.csv.gz',
    '/work/shah/tantalus/SC-1561/results/results/alignment/A73044A_alignment_metrics.csv.gz',
    '/work/shah/tantalus/SC-1583/results/results/alignment/A96139A_alignment_metrics.csv.gz',
    '/work/shah/tantalus/SC-1585/results/results/alignment/A96147A_alignment_metrics.csv.gz',
]

rule all_fig4:
    input:
        expand(
            'plots/fig4/{dataset}/cn_heatmaps.png',
            dataset=[d for d in datasets]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[d for d in datasets]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro_filtered.png',
            dataset=[d for d in datasets]
        ),
        'plots/fig4/all/rt_corr.png',


# fetch the raw data
rule get_data_4:
    input:
        cn_data_urls = cn_data_urls,
        metrics_data_urls = metrics_data_urls,
        align_metrics_data_urls = align_metrics_data_urls,
    output:
        cn_data = 'analysis/fig4/all/cn_data.tsv',
        metrics_data = 'analysis/fig4/all/metrics_data.tsv',
        align_metrics_data = 'analysis/fig4/all/align_metrics_data.tsv'
    log: 'logs/fig4/all/get_data.log'
    shell:
        'python scripts/fig4/get_data.py '
        '--cn_data_urls {input.cn_data_urls} '
        '--metrics_data_urls {input.metrics_data_urls} '
        '--align_metrics_data_urls {input.align_metrics_data_urls} '
        '--cn_out {output.cn_data} '
        '--metrics_out {output.metrics_data} '
        '--align_metrics_out {output.align_metrics_data} '
        '&> {log}'


# make sure all cells have same loci and no NaNs
rule filter_data_4:
    input: 
        cn_input = 'analysis/fig4/all/cn_data.tsv',
        metrics_data = 'analysis/fig4/all/metrics_data.tsv'
    output: 'analysis/fig4/all/cn_data_filtered.tsv'
    log: 'logs/fig4/all/filter_data.log'
    shell:
        'python scripts/fig4/filter_data.py '
        '{input} {params} {output} &> {log}'


rule compute_ccc_features_4:
    input: 'analysis/fig4/all/cn_data_filtered.tsv'
    output: 'analysis/fig4/all/cn_data_features.tsv'
    log: 'logs/fig4/all/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# split the cn and metrics data by sample_id which corresponds to cell line
rule split_cell_line_4:
    input:
        cn_data = 'analysis/fig4/all/cn_data_features.tsv',
        metrics_data = 'analysis/fig4/all/metrics_data.tsv'
    output:
        cn_T47D = 'analysis/fig4/T47D/cn_data_features.tsv',
        metrics_T47D = 'analysis/fig4/T47D/metrics_data.tsv',
        cn_GM18507 = 'analysis/fig4/GM18507/cn_data_features.tsv',
        metrics_GM18507 = 'analysis/fig4/GM18507/metrics_data.tsv',
    run:
        cn = pd.read_csv(str(input.cn_data), sep='\t')
        metrics = pd.read_csv(str(input.metrics_data), sep='\t')

        cn_T47D = cn.query("sample_id=='SA1044'")
        metrics_T47D = metrics.query("sample_id=='SA1044'")
        cn_GM18507 = cn.query("sample_id=='SA928'")
        metrics_GM18507 = metrics.query("sample_id=='SA928'")

        cn_T47D.to_csv(str(output.cn_T47D), sep='\t', index=False)
        metrics_T47D.to_csv(str(output.metrics_T47D), sep='\t', index=False)
        cn_GM18507.to_csv(str(output.cn_GM18507), sep='\t', index=False)
        metrics_GM18507.to_csv(str(output.metrics_GM18507), sep='\t', index=False)


# use metrics file to split each cell in filtered cn data by cell cycle state
rule split_cell_cycle_4:
    input: 'analysis/fig4/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/fig4/{dataset}/cn_s.tsv',
        cn_g1 = 'analysis/fig4/{dataset}/cn_g1.tsv',
        cn_g2 = 'analysis/fig4/{dataset}/cn_g2.tsv'
    log: 'logs/fig4/{dataset}/split_cell_cycle.log'
    shell:
        'python scripts/fig4/split_cell_cycle.py '
        '{input} {params} {output} &> {log}'


rule plot_cn_heatmaps_4:
    input:
        s_phase = 'analysis/fig4/{dataset}/cn_s.tsv',
        g1_phase = 'analysis/fig4/{dataset}/cn_g1.tsv',
        g2_phase = 'analysis/fig4/{dataset}/cn_g2.tsv',
    output: 'plots/fig4/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = 'flow-sorted'
    log:
        'logs/fig4/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig4/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule infer_scRT_pyro_4:
    input:
        cn_s = 'analysis/fig4/{dataset}/cn_s.tsv',
        cn_g1 = 'analysis/fig4/{dataset}/cn_g1.tsv',
        cn_g2 = 'analysis/fig4/{dataset}/cn_g2.tsv'
    output: 'analysis/fig4/{dataset}/cn_s_pyro_infered.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_cells',
        infer_mode = 'pyro'
    log: 'logs/fig4/{dataset}/infer_scRT_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_4:
    input: 'analysis/fig4/{dataset}/cn_s_pyro_infered.tsv'
    output: 
        plot1 = 'plots/fig4/{dataset}/scRT_heatmaps_pyro.png',
        plot2 = 'plots/fig4/{dataset}/frac_rt_distributions_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'model_s_time'
    log: 'logs/fig4/{dataset}/plot_inferred_cn_vs_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule remove_nonreplicating_cells_4:
    input: 
        cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered.tsv',
        cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered.tsv',
        cn_all = 'analysis/fig4/all/cn_s_pyro_infered.tsv'
    output:
        cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered_filtered.tsv',
        cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered_filtered.tsv',
        cn_all = 'analysis/fig4/all/cn_s_pyro_infered_filtered.tsv'
    params:
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/fig4/remove_nonreplicating_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/remove_nonreplicating_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_filtered_4:
    input: 'analysis/fig4/{dataset}/cn_s_pyro_infered_filtered.tsv'
    output: 
        plot1 = 'plots/fig4/{dataset}/scRT_heatmaps_pyro_filtered.png',
        plot2 = 'plots/fig4/{dataset}/frac_rt_distributions_pyro_filtered.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/fig4/{dataset}/plot_inferred_cn_vs_scRT_filtered.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_4:
    input: 
        cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered_filtered.tsv',
        cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered_filtered.tsv',
        cn_all = 'analysis/fig4/all/cn_s_pyro_infered_filtered.tsv'
    output: 'analysis/fig4/all/rt_pseudobulks.tsv'
    log: 'logs/fig4/all/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_rt_profiles_4:
    input: 'analysis/fig4/all/rt_pseudobulks.tsv'
    output:
        plot1 = 'plots/fig4/all/rt_diff_split.png',
        plot2 = 'plots/fig4/all/rt_diff_joint.png',
        plot3 = 'plots/fig4/all/rt_corr.png',
        plot4 = 'plots/fig4/all/rt_split_chr1.png',
        plot5 = 'plots/fig4/all/rt_joint_chr1.png',
    log: 'logs/fig4/all/plot_rt_profiles.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/plot_rt_profiles.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# TODO: update this rule for computing T-width for the joint vs split versions of the model
# rule twidth_analysis_pyro_4:
#     input: 
#         'analysis/fig2/{dataset}/s_phase_cells_pyro_infered.tsv'
#     output: 
#         plot1 = 'plots/fig2/{dataset}/twidth_heatmaps_pyro.png',
#         plot2 = 'plots/fig2/{dataset}/twidth_curves_pyro.png',
#     params:
#         dataset = lambda wildcards: wildcards.dataset,
#         A = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['A'],
#         nb_r = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['nb_r'],
#         rt_col = lambda wildcards: config['simulated_datasets'][wildcards.dataset]['rt_col'],
#         frac_rt_col = 'model_s_time',
#         true_frac_col = 'true_t',
#         rep_state = 'model_rep_state',
#         true_rep_state = 'true_rep',
#         infer_mode = 'pyro'
#     log: 'logs/fig2/{dataset}/twidth_analysis_pyro.log'
#     shell:
#         'source ../scgenome/venv/bin/activate ; '
#         'python3 scripts/fig2/twidth_analysis.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'

