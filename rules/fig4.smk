import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

datasets = ['all', 'GM18507', 'T47D']

# create list of datasets that also includes the permuted datasets
perm_datasets = [y for x in [datasets, config['permuted_datasets']] for y in x]

cn_data_urls = [
    '/juno/work/shah/isabl_data_lake/experiments/78/11/17811/analyses/scdna-hmmcopy__0.0.4__27278/results/A90553C_reads.csv.gz',
    '/juno/work/shah/isabl_data_lake/experiments/71/34/17134/analyses/scdna-hmmcopy__0.0.3__22933/results/A73044A_reads.csv.gz',
    '/juno/work/shah/isabl_data_lake/experiments/78/12/17812/analyses/scdna-hmmcopy__0.0.4__27279/results/SA1044_reads.csv.gz',  # A96139A
    '/juno/work/shah/isabl_data_lake/experiments/78/13/17813/analyses/scdna-hmmcopy__0.0.4__27280/results/SA1044_reads.csv.gz',  # A96147A
]

metrics_data_urls = [
    '/juno/work/shah/isabl_data_lake/experiments/78/11/17811/analyses/scdna-annotation__0.0.4__27286/results/A90553C_metrics.csv.gz',
    '/juno/work/shah/isabl_data_lake/experiments/71/34/17134/analyses/scdna-annotation__0.0.3__22937/results/A73044A_metrics.csv.gz',
    '/juno/work/shah/isabl_data_lake/experiments/78/12/17812/analyses/scdna-annotation__0.0.4__27287/results/SA1044_metrics.csv.gz',  # A96139A
    '/juno/work/shah/isabl_data_lake/experiments/78/13/17813/analyses/scdna-annotation__0.0.4__27288/results/SA1044_metrics.csv.gz',  # A96147A
]

align_metrics_data_urls = [
    '/juno/work/shah/isabl_data_lake/experiments/78/11/17811/analyses/scdna-alignment__0.0.4__27195/results/A90553C_alignment_metrics.csv.gz',
    '/juno/work/shah/isabl_data_lake/experiments/71/34/17134/analyses/scdna-alignment__0.0.3__22930/results/A73044A_alignment_metrics.csv.gz',
    '/juno/work/shah/isabl_data_lake/experiments/78/12/17812/analyses/scdna-alignment__0.0.4__27196/results/SA1044_alignment_metrics.csv.gz',  # A96139A
    '/juno/work/shah/isabl_data_lake/experiments/78/13/17813/analyses/scdna-alignment__0.0.4__27197/results/SA1044_alignment_metrics.csv.gz',  # A96147A
]

rule all_fig4:
    input:
        expand(
            'plots/fig4/{dataset}/cn_heatmaps.png',
            dataset=[d for d in datasets]
        ),
        expand(
            'plots/fig4/{dataset}/cn_clone_heatmaps.png',
            dataset=[d for d in perm_datasets]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[d for d in datasets]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro_composite.png',
            dataset=[d for d in perm_datasets]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro_filtered.png',
            dataset=[d for d in datasets]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro_composite_filtered.png',
            dataset=[d for d in perm_datasets]
        ),
        expand(
            'analysis/fig4/{dataset}/rt_pseudobulks_composite.tsv',
            dataset=[
                d for d in perm_datasets
                if (d not in ['T47D', 'GM18507'])
            ]
        ),
        'plots/fig4/all/rt_corr.png',
        'plots/fig4/all/rt_corr_composite.png',
        'plots/fig4/all/twidth_curves.png',
        'plots/fig4/all/twidth_curves_composite.png',
        'plots/fig4/permuted/summary.png',
        'plots/fig4/permuted/rt_corr_composite.png'


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


# swaps the cell_cycle_state for x% of the G1/2-phase cells and pass into model
# these datasets should evaluate how well the model performs as a classifier
rule permute_cell_cycle_labels_4:
    input: 'analysis/fig4/all/cn_data_features.tsv'
    output: 'analysis/fig4/{dataset}/cn_data_features.tsv'
    params:
        permute_rate = lambda wildcards: config['permuted_datasets'][wildcards.dataset]['rate'],
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig4/{dataset}/permute_cell_cycle_labels.log'
    shell:
        'python scripts/fig4/permute_cell_cycle_labels.py '
        '{input} {params} {output} &> {log}'


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


# cluster all the g1/2-phase cells into clones (using K-means instead of sitka for now)
rule cluster_into_clones_4:
    input:
        cn_g1 = 'analysis/fig4/{dataset}/cn_g1.tsv',
        cn_g2 = 'analysis/fig4/{dataset}/cn_g2.tsv'
    output: 'analysis/fig4/{dataset}/cn_g_with_clone_id.tsv'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig4/{dataset}/cluster_into_clones.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python scripts/fig4/cluster_into_clones.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_cn_clone_heatmaps_4:
    input:
        s_phase = 'analysis/fig4/{dataset}/cn_s.tsv',
        g_phase = 'analysis/fig4/{dataset}/cn_g_with_clone_id.tsv',
    output: 'plots/fig4/{dataset}/cn_clone_heatmaps.png'
    params:
        value_col = 'state',
        dataset = 'flow-sorted'
    log:
        'logs/fig4/{dataset}/plot_cn_clone_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig4/plot_cn_clone_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule infer_scRT_pyro_4:
    input:
        cn_s = 'analysis/fig4/{dataset}/cn_s.tsv',
        cn_g = 'analysis/fig4/{dataset}/cn_g_with_clone_id.tsv'
    output: 'analysis/fig4/{dataset}/cn_s_pyro_infered.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_clones',
    log: 'logs/fig4/{dataset}/infer_scRT_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_composite_4:
    input:
        cn_s = 'analysis/fig4/{dataset}/cn_s.tsv',
        cn_g = 'analysis/fig4/{dataset}/cn_g_with_clone_id.tsv'
    output: 'analysis/fig4/{dataset}/cn_s_pyro_composite_infered.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
    log: 'logs/fig4/{dataset}/infer_scRT_pyro_composite.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/infer_scRT_composite.py '
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


rule plot_inferred_cn_vs_scRT_composite_4:
    input: 'analysis/fig4/{dataset}/cn_s_pyro_composite_infered.tsv'
    output: 
        plot1 = 'plots/fig4/{dataset}/scRT_heatmaps_pyro_composite.png',
        plot2 = 'plots/fig4/{dataset}/frac_rt_distributions_pyro_composite.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'model_s_time'
    log: 'logs/fig4/{dataset}/plot_inferred_cn_vs_scRT_composite.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# rule remove_nonreplicating_cells_4:
#     input: 
#         cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered.tsv',
#         cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered.tsv',
#         cn_all = 'analysis/fig4/all/cn_s_pyro_infered.tsv'
#     output:
#         cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered_filtered.tsv',
#         cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered_filtered.tsv',
#         cn_all = 'analysis/fig4/all/cn_s_pyro_infered_filtered.tsv'
#     params:
#         frac_rt_col = 'cell_frac_rep'
#     log: 'logs/fig4/remove_nonreplicating_cells.log'
#     shell:
#         'source ../scdna_replication_tools/venv/bin/activate ; '
#         'python3 scripts/fig4/remove_nonreplicating_cells.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'


rule remove_nonreplicating_cells_simple_4:
    input: 'analysis/fig4/{dataset}/cn_s_pyro_infered.tsv',
    output:
        good = 'analysis/fig4/{dataset}/cn_s_pyro_infered_filtered.tsv',
        bad = 'analysis/fig4/{dataset}/model_nonrep_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state'
    log: 'logs/fig4/{dataset}/remove_nonreplicating_cells_simple.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/remove_nonreplicating_cells_simple.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# rule remove_nonreplicating_cells_composite_4:
#     input: 
#         cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_composite_infered.tsv',
#         cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_composite_infered.tsv',
#         cn_all = 'analysis/fig4/all/cn_s_pyro_composite_infered.tsv'
#     output:
#         cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered_composite_filtered.tsv',
#         cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered_composite_filtered.tsv',
#         cn_all = 'analysis/fig4/all/cn_s_pyro_infered_composite_filtered.tsv'
#     params:
#         frac_rt_col = 'cell_frac_rep'
#     log: 'logs/fig4/remove_nonreplicating_cells_composite.log'
#     shell:
#         'source ../scdna_replication_tools/venv/bin/activate ; '
#         'python3 scripts/fig4/remove_nonreplicating_cells.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'


rule remove_nonreplicating_cells_composite_simple_4:
    input: 'analysis/fig4/{dataset}/cn_s_pyro_composite_infered.tsv',
    output:
        good = 'analysis/fig4/{dataset}/cn_s_pyro_infered_composite_filtered.tsv',
        bad = 'analysis/fig4/{dataset}/model_nonrep_cells_composite.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state'
    log: 'logs/fig4/{dataset}/remove_nonreplicating_cells_composite_simple.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/remove_nonreplicating_cells_simple.py '
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


rule plot_inferred_cn_vs_scRT_composite_filtered_4:
    input: 'analysis/fig4/{dataset}/cn_s_pyro_infered_composite_filtered.tsv'
    output: 
        plot1 = 'plots/fig4/{dataset}/scRT_heatmaps_pyro_composite_filtered.png',
        plot2 = 'plots/fig4/{dataset}/frac_rt_distributions_pyro_composite_filtered.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/fig4/{dataset}/plot_inferred_cn_vs_scRT_composite_filtered.log'
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


rule compute_rt_pseudobulks_composite_4:
    input: 
        cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered_composite_filtered.tsv',
        cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered_composite_filtered.tsv',
        cn_all = 'analysis/fig4/all/cn_s_pyro_infered_composite_filtered.tsv'
    output: 'analysis/fig4/all/rt_pseudobulks_composite.tsv'
    log: 'logs/fig4/all/compute_rt_pseudobulks_composite.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_permuted_composite_4:
    input: 'analysis/fig4/{dataset}/cn_s_pyro_infered_composite_filtered.tsv',
    output: 'analysis/fig4/{dataset}/rt_pseudobulks_composite.tsv'
    params:
        rep_col = 'model_rep_state'
    log: 'logs/fig4/{dataset}/compute_rt_pseudobulks_composite.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/compute_rt_pseudobulks_simple.py '
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


rule plot_rt_profiles_composite_4:
    input: 'analysis/fig4/all/rt_pseudobulks_composite.tsv'
    output:
        plot1 = 'plots/fig4/all/rt_diff_split_composite.png',
        plot2 = 'plots/fig4/all/rt_diff_joint_composite.png',
        plot3 = 'plots/fig4/all/rt_corr_composite.png',
        plot4 = 'plots/fig4/all/rt_split_chr1_composite.png',
        plot5 = 'plots/fig4/all/rt_joint_chr1_composite.png',
    log: 'logs/fig4/all/plot_rt_profiles_composite.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/plot_rt_profiles.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# def get_permuted_cn_good(wildcards):
#     files = []
#     for d in config['permuted_datasets']:
#         files.append('analysis/fig4/{}/cn_s_pyro_infered_composite_filtered.tsv'.format(d))
#     return expand(files)


# def get_permuted_cn_bad(wildcards):
#     files = []
#     for d in config['permuted_datasets']:
#         files.append('analysis/fig4/{}/model_nonrep_cells_composite.tsv'.format(d))
#     return expand(files)


rule analyze_permuted_datasets_4:
    input: 
        cn_good = expand(
            'analysis/fig4/{dataset}/cn_s_pyro_infered_composite_filtered.tsv',
            dataset=[d for d in config['permuted_datasets']]
        ),
        cn_bad = expand(
            'analysis/fig4/{dataset}/model_nonrep_cells_composite.tsv',
            dataset=[d for d in config['permuted_datasets']]
        )
    output:
        summary = 'analysis/fig4/permuted/summary.tsv',
        cell_metrics = 'analysis/fig4/permuted/cell_metrics.tsv',
        summary_plots = 'plots/fig4/permuted/summary.png',
        ccc_plots = 'plots/fig4/permuted/ccc_features.png',
    params:
        datasets = expand([d for d in config['permuted_datasets']]),
        rates = expand([str(config['permuted_datasets'][d]['rate']) for d in config['permuted_datasets']])
    log: 'logs/fig4/permuted/analyze_permuted_datasets.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/analyze_permuted_datasets.py '
        '--cn_good {input.cn_good} '
        '--cn_bad {input.cn_bad} '
        '--datasets {params.datasets} '
        '--rates {params.rates} '
        '--summary_output {output.summary} '
        '--metrics_output {output.cell_metrics} '
        '--summary_plots {output.summary_plots} '
        '--ccc_plots {output.ccc_plots} '
        '&> {log} ; '
        'deactivate'


rule permuted_dataset_rt_profiles_4:
    input: 
        rt_all = 'analysis/fig4/all/rt_pseudobulks_composite.tsv',
        rt_perm = expand(
            'analysis/fig4/{dataset}/model_nonrep_cells_composite.tsv',
            dataset=[d for d in config['permuted_datasets']]
        )
    output:
        rt_table = 'analysis/fig4/permuted/rt_pseudobulks_composite.tsv',
        cor_plot = 'plots/fig4/permuted/rt_corr_composite.png',
    params:
        datasets = expand([d for d in config['permuted_datasets']])
    log: 'logs/fig4/permuted/permuted_dataset_rt_profiles.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/permuted_dataset_rt_profiles.py '
        '--rt_all {input.rt_all} '
        '--rt_perm {input.rt_perm} '
        '--datasets {params.datasets} '
        '--rt_table {output.rt_table} '
        '--cor_plot {output.cor_plot} '
        '&> {log} ; '
        'deactivate'


rule twidth_analysis_pyro_4:
    input: 
        cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered_filtered.tsv',
        cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered_filtered.tsv',
        cn_all = 'analysis/fig4/all/cn_s_pyro_infered_filtered.tsv',
        pseudobulk = 'analysis/fig4/all/rt_pseudobulks.tsv'
    output: 
        output_tsv = 'analysis/fig4/all/twidth_values.tsv',
        output_png = 'plots/fig4/all/twidth_curves.png',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_state = 'model_rep_state'
    log: 'logs/fig4/all/twidth_analysis_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_pyro_composite_4:
    input: 
        cn_T47D = 'analysis/fig4/T47D/cn_s_pyro_infered_composite_filtered.tsv',
        cn_GM18507 = 'analysis/fig4/GM18507/cn_s_pyro_infered_composite_filtered.tsv',
        cn_all = 'analysis/fig4/all/cn_s_pyro_infered_composite_filtered.tsv',
        pseudobulk = 'analysis/fig4/all/rt_pseudobulks_composite.tsv'
    output: 
        output_tsv = 'analysis/fig4/all/twidth_values_composite.tsv',
        output_png = 'plots/fig4/all/twidth_curves_composite.png',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_state = 'model_rep_state'
    log: 'logs/fig4/all/twidth_analysis_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'

