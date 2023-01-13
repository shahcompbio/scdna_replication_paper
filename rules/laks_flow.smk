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

rule all_laks_flow:
    input:
        expand(
            'plots/laks_flow/{dataset}/cn_heatmaps.png',
            dataset=[d for d in datasets]
        ),
        expand(
            'plots/laks_flow/{dataset}/cn_clone_heatmaps.png',
            dataset=[d for d in perm_datasets]
        ),
        # expand(
        #     'plots/laks_flow/{dataset}/scRT_heatmaps_pyro.png',
        #     dataset=[d for d in datasets]
        # ),
        expand(
            'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_composite.png',
            dataset=[d for d in perm_datasets]
        ),
        # expand(
        #     'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_g.png',
        #     dataset=[d for d in perm_datasets]
        # ),
        # expand(
        #     'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_filtered.png',
        #     dataset=[d for d in datasets]
        # ),
        expand(
            'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_composite_filtered.png',
            dataset=[d for d in perm_datasets]
        ),
        # expand(
        #     'analysis/laks_flow/{dataset}/rt_pseudobulks_composite.tsv',
        #     dataset=[d for d in ['T47D']]
        # ),
        expand(
            'analysis/laks_flow/{dataset}/rt_pseudobulks_composite.tsv',
            dataset=[
                d for d in perm_datasets
                if (d not in ['T47D', 'GM18507'])
            ]
        ),
        expand(
            'analysis/laks_flow/{dataset}/cn_pseudobulks.tsv',
            dataset=[
                d for d in ['T47D', 'GM18507', 'all']
            ]
        ),
        'plots/laks_flow/GM18507/cn_s_example.png',
        'plots/laks_flow/all/flow_error_cells.png',
        'plots/laks_flow/all/rpm_umap.png',
        'plots/laks_flow/all/rt_corr.png',
        'plots/laks_flow/all/rt_corr_composite.png',
        'plots/laks_flow/all/twidth_curves.png',
        'plots/laks_flow/all/twidth_curves_composite.png',
        'plots/laks_flow/permuted/ccc_features.png',
        'plots/laks_flow/permuted/summary.png',
        'plots/laks_flow/permuted/rt_corr_composite.png'


# fetch the raw data
rule get_data_lf:
    input:
        cn_data_urls = cn_data_urls,
        metrics_data_urls = metrics_data_urls,
        align_metrics_data_urls = align_metrics_data_urls,
    output:
        cn_data = 'analysis/laks_flow/all/cn_data.tsv',
        metrics_data = 'analysis/laks_flow/all/metrics_data.tsv',
        align_metrics_data = 'analysis/laks_flow/all/align_metrics_data.tsv'
    log: 'logs/laks_flow/all/get_data.log'
    shell:
        'python scripts/laks_flow/get_data.py '
        '--cn_data_urls {input.cn_data_urls} '
        '--metrics_data_urls {input.metrics_data_urls} '
        '--align_metrics_data_urls {input.align_metrics_data_urls} '
        '--cn_out {output.cn_data} '
        '--metrics_out {output.metrics_data} '
        '--align_metrics_out {output.align_metrics_data} '
        '&> {log}'


# make sure all cells have same loci and no NaNs
rule filter_data_lf:
    input: 
        cn_input = 'analysis/laks_flow/all/cn_data.tsv',
        metrics_data = 'analysis/laks_flow/all/metrics_data.tsv'
    output: 'analysis/laks_flow/all/cn_data_filtered.tsv'
    log: 'logs/laks_flow/all/filter_data.log'
    shell:
        'python scripts/laks_flow/filter_data.py '
        '{input} {params} {output} &> {log}'


rule compute_ccc_features_lf:
    input: 'analysis/laks_flow/all/cn_data_filtered.tsv'
    output: 'analysis/laks_flow/all/cn_data_features.tsv'
    log: 'logs/laks_flow/all/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# split the cn and metrics data by sample_id which corresponds to cell line
rule split_cell_line_lf:
    input:
        cn_data = 'analysis/laks_flow/all/cn_data_features.tsv',
        metrics_data = 'analysis/laks_flow/all/metrics_data.tsv'
    output:
        cn_T47D = 'analysis/laks_flow/T47D/cn_data_features.tsv',
        metrics_T47D = 'analysis/laks_flow/T47D/metrics_data.tsv',
        cn_GM18507 = 'analysis/laks_flow/GM18507/cn_data_features.tsv',
        metrics_GM18507 = 'analysis/laks_flow/GM18507/metrics_data.tsv',
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
rule permute_cell_cycle_labels_lf:
    input: 'analysis/laks_flow/all/cn_data_features.tsv'
    output: 'analysis/laks_flow/{dataset}/cn_data_features.tsv'
    params:
        permute_rate = lambda wildcards: config['permuted_datasets'][wildcards.dataset]['rate'],
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/laks_flow/{dataset}/permute_cell_cycle_labels.log'
    shell:
        'python scripts/laks_flow/permute_cell_cycle_labels.py '
        '{input} {params} {output} &> {log}'


# use metrics file to split each cell in filtered cn data by cell cycle state
rule split_cell_cycle_lf:
    input: 'analysis/laks_flow/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/laks_flow/{dataset}/cn_s.tsv',
        cn_g1 = 'analysis/laks_flow/{dataset}/cn_g1.tsv',
        cn_g2 = 'analysis/laks_flow/{dataset}/cn_g2.tsv'
    log: 'logs/laks_flow/{dataset}/split_cell_cycle.log'
    shell:
        'python scripts/laks_flow/split_cell_cycle.py '
        '{input} {params} {output} &> {log}'


rule plot_cn_heatmaps_lf:
    input:
        s_phase = 'analysis/laks_flow/{dataset}/cn_s.tsv',
        g1_phase = 'analysis/laks_flow/{dataset}/cn_g1.tsv',
        g2_phase = 'analysis/laks_flow/{dataset}/cn_g2.tsv',
    output: 'plots/laks_flow/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = 'flow-sorted'
    log:
        'logs/laks_flow/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/laks_flow/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


# cluster all the g1/2-phase cells into clones (using K-means instead of sitka for now)
rule cluster_into_clones_lf:
    input:
        cn_g1 = 'analysis/laks_flow/{dataset}/cn_g1.tsv',
        cn_g2 = 'analysis/laks_flow/{dataset}/cn_g2.tsv'
    output: 'analysis/laks_flow/{dataset}/cn_g_with_clone_id.tsv'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/laks_flow/{dataset}/cluster_into_clones.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python scripts/laks_flow/cluster_into_clones.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_cn_clone_heatmaps_lf:
    input:
        s_phase = 'analysis/laks_flow/{dataset}/cn_s.tsv',
        g_phase = 'analysis/laks_flow/{dataset}/cn_g_with_clone_id.tsv',
    output: 'plots/laks_flow/{dataset}/cn_clone_heatmaps.png'
    params:
        value_col = 'state',
        dataset = 'flow-sorted'
    log:
        'logs/laks_flow/{dataset}/plot_cn_clone_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/laks_flow/plot_cn_clone_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule infer_scRT_pyro_lf:
    input:
        cn_s = 'analysis/laks_flow/{dataset}/cn_s.tsv',
        cn_g = 'analysis/laks_flow/{dataset}/cn_g_with_clone_id.tsv'
    output:
        main_s_out = 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred.tsv',
        supp_s_out = 'analysis/laks_flow/{dataset}/scRT_pyro_supp_s_output.tsv',
        main_g_out = 'analysis/laks_flow/{dataset}/cn_g_pyro_inferred.tsv',
        supp_g_out = 'analysis/laks_flow/{dataset}/scRT_pyro_supp_g_output.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_clones',
    log: 'logs/laks_flow/{dataset}/infer_scRT_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_composite_lf:
    input:
        cn_s = 'analysis/laks_flow/{dataset}/cn_s.tsv',
        cn_g = 'analysis/laks_flow/{dataset}/cn_g_with_clone_id.tsv'
    output:
        main_s_out = 'analysis/laks_flow/{dataset}/cn_s_pyro_composite_inferred.tsv',
        supp_s_out = 'analysis/laks_flow/{dataset}/scRT_pyro_composite_supp_s_output.tsv',
        main_g_out = 'analysis/laks_flow/{dataset}/cn_g_pyro_composite_inferred.tsv',
        supp_g_out = 'analysis/laks_flow/{dataset}/scRT_pyro_composite_supp_g_output.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
    log: 'logs/laks_flow/{dataset}/infer_scRT_pyro_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/infer_scRT_composite.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_lf:
    input: 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred.tsv'
    output: 
        plot1 = 'plots/laks_flow/{dataset}/scRT_heatmaps_pyro.png',
        plot2 = 'plots/laks_flow/{dataset}/frac_rt_distributions_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'model_tau'
    log: 'logs/laks_flow/{dataset}/plot_inferred_cn_vs_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_g_vs_scRT_lf:
    input: 'analysis/laks_flow/{dataset}/cn_g_pyro_inferred.tsv'
    output: 
        plot1 = 'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_g.png',
        plot2 = 'plots/laks_flow/{dataset}/frac_rt_distributions_pyro_g.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'model_tau'
    log: 'logs/laks_flow/{dataset}/plot_inferred_cn_g_vs_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_composite_lf:
    input: 'analysis/laks_flow/{dataset}/cn_s_pyro_composite_inferred.tsv'
    output: 
        plot1 = 'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_composite.png',
        plot2 = 'plots/laks_flow/{dataset}/frac_rt_distributions_pyro_composite.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'model_tau'
    log: 'logs/laks_flow/{dataset}/plot_inferred_cn_vs_scRT_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_g_vs_scRT_composite_lf:
    input: 'analysis/laks_flow/{dataset}/cn_g_pyro_composite_inferred.tsv'
    output: 
        plot1 = 'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_composite_g.png',
        plot2 = 'plots/laks_flow/{dataset}/frac_rt_distributions_pyro_composite_g.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'model_tau'
    log: 'logs/laks_flow/{dataset}/plot_inferred_cn_g_vs_scRT_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule revise_cell_cycle_labels_lf:
    input: 
        cn_s = 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred.tsv',
        cn_g = 'analysis/laks_flow/{dataset}/cn_g_pyro_inferred.tsv',
    output:
        out_s = 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred_filtered.tsv',
        out_g = 'analysis/laks_flow/{dataset}/cn_g_pyro_inferred_filtered.tsv',
        out_lowqual = 'analysis/laks_flow/{dataset}/cn_lowqual_pyro_inferred_filtered.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/laks_flow/{dataset}/revise_cell_cycle_labels.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/revise_cell_cycle_labels.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule revise_cell_cycle_labels_composite_lf:
    input: 
        cn_s = 'analysis/laks_flow/{dataset}/cn_s_pyro_composite_inferred.tsv',
        cn_g = 'analysis/laks_flow/{dataset}/cn_g_pyro_composite_inferred.tsv',
    output:
        out_s = 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred_composite_filtered.tsv',
        out_g = 'analysis/laks_flow/{dataset}/cn_g_pyro_inferred_composite_filtered.tsv',
        out_lowqual = 'analysis/laks_flow/{dataset}/cn_lowqual_pyro_inferred_composite_filtered.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/laks_flow/{dataset}/revise_cell_cycle_labels_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/revise_cell_cycle_labels.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# TODO: figure out why I can't load the data for this rule
rule rpm_umap_lf:
    input: 
        cn_s = 'analysis/laks_flow/all/cn_s_pyro_inferred_composite_filtered.tsv',
        cn_g = 'analysis/laks_flow/all/cn_g_pyro_inferred_composite_filtered.tsv',
    output: 'plots/laks_flow/all/rpm_umap.png'
    params:
        rpm_col = 'rpm'
    log: 'logs/laks_flow/all/rpm_umap.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/rpm_umap.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_filtered_lf:
    input: 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred_filtered.tsv'
    output: 
        plot1 = 'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_filtered.png',
        plot2 = 'plots/laks_flow/{dataset}/frac_rt_distributions_pyro_filtered.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/laks_flow/{dataset}/plot_inferred_cn_vs_scRT_filtered.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_composite_filtered_lf:
    input: 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred_composite_filtered.tsv'
    output: 
        plot1 = 'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_composite_filtered.png',
        plot2 = 'plots/laks_flow/{dataset}/frac_rt_distributions_pyro_composite_filtered.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/laks_flow/{dataset}/plot_inferred_cn_vs_scRT_composite_filtered.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_cn_pseudobulks_lf:
    input: 'analysis/laks_flow/{dataset}/cn_g_with_clone_id.tsv'
    output: 'analysis/laks_flow/{dataset}/cn_pseudobulks.tsv'
    params:
        cn_state_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/laks_flow/{dataset}/compute_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/compute_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_lf:
    input: 
        cn_T47D = 'analysis/laks_flow/T47D/cn_s_pyro_inferred_filtered.tsv',
        cn_GM18507 = 'analysis/laks_flow/GM18507/cn_s_pyro_inferred_filtered.tsv',
        cn_all = 'analysis/laks_flow/all/cn_s_pyro_inferred_filtered.tsv'
    output: 'analysis/laks_flow/all/rt_pseudobulks.tsv'
    log: 'logs/laks_flow/all/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_composite_lf:
    input: 
        cn_T47D = 'analysis/laks_flow/T47D/cn_s_pyro_inferred_composite_filtered.tsv',
        cn_GM18507 = 'analysis/laks_flow/GM18507/cn_s_pyro_inferred_composite_filtered.tsv',
        cn_all = 'analysis/laks_flow/all/cn_s_pyro_inferred_composite_filtered.tsv'
    output: 'analysis/laks_flow/all/rt_pseudobulks_composite.tsv'
    log: 'logs/laks_flow/all/compute_rt_pseudobulks_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_permuted_composite_lf:
    input: 'analysis/laks_flow/{dataset}/cn_s_pyro_inferred_composite_filtered.tsv',
    output: 'analysis/laks_flow/{dataset}/rt_pseudobulks_composite.tsv'
    params:
        rep_col = 'model_rep_state'
    log: 'logs/laks_flow/{dataset}/compute_rt_pseudobulks_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/compute_rt_pseudobulks_simple.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_rt_profiles_lf:
    input: 'analysis/laks_flow/all/rt_pseudobulks.tsv'
    output:
        plot1 = 'plots/laks_flow/all/rt_diff_split.png',
        plot2 = 'plots/laks_flow/all/rt_diff_joint.png',
        plot3 = 'plots/laks_flow/all/rt_corr.png',
        plot4 = 'plots/laks_flow/all/rt_split_chr1.png',
        plot5 = 'plots/laks_flow/all/rt_joint_chr1.png',
    log: 'logs/laks_flow/all/plot_rt_profiles.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_rt_profiles.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_rt_profiles_composite_lf:
    input: 'analysis/laks_flow/all/rt_pseudobulks_composite.tsv'
    output:
        plot1 = 'plots/laks_flow/all/rt_diff_split_composite.png',
        plot2 = 'plots/laks_flow/all/rt_diff_joint_composite.png',
        plot3 = 'plots/laks_flow/all/rt_corr_composite.png',
        plot4 = 'plots/laks_flow/all/rt_split_chr1_composite.png',
        plot5 = 'plots/laks_flow/all/rt_joint_chr1_composite.png',
    log: 'logs/laks_flow/all/plot_rt_profiles_composite.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_rt_profiles.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# def get_permuted_cn_good(wildcards):
#     files = []
#     for d in config['permuted_datasets']:
#         files.append('analysis/laks_flow/{}/cn_s_pyro_inferred_composite_filtered.tsv'.format(d))
#     return expand(files)


# def get_permuted_cn_bad(wildcards):
#     files = []
#     for d in config['permuted_datasets']:
#         files.append('analysis/laks_flow/{}/model_nonrep_cells_composite.tsv'.format(d))
#     return expand(files)


rule analyze_permuted_datasets_lf:
    input: 
        cn_s = expand(
            'analysis/laks_flow/{dataset}/cn_s_pyro_inferred_composite_filtered.tsv',
            dataset=[d for d in config['permuted_datasets']]
        ),
        cn_g = expand(
            'analysis/laks_flow/{dataset}/cn_g_pyro_inferred_composite_filtered.tsv',
            dataset=[d for d in config['permuted_datasets']]
        )
    output:
        summary = 'analysis/laks_flow/permuted/summary.tsv',
        cell_metrics = 'analysis/laks_flow/permuted/cell_metrics.tsv',
        summary_plots = 'plots/laks_flow/permuted/summary.png',
        ccc_plots = 'plots/laks_flow/permuted/ccc_features.png',
    params:
        datasets = expand([d for d in config['permuted_datasets']]),
        rates = expand([str(config['permuted_datasets'][d]['rate']) for d in config['permuted_datasets']])
    log: 'logs/laks_flow/permuted/analyze_permuted_datasets.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/analyze_permuted_datasets.py '
        '--cn_good {input.cn_s} '
        '--cn_bad {input.cn_g} '
        '--datasets {params.datasets} '
        '--rates {params.rates} '
        '--summary_output {output.summary} '
        '--metrics_output {output.cell_metrics} '
        '--summary_plots {output.summary_plots} '
        '--ccc_plots {output.ccc_plots} '
        '&> {log} ; '
        'deactivate'


rule permuted_dataset_rt_profiles_lf:
    input: 
        rt_ref = 'analysis/laks_flow/all/rt_pseudobulks_composite.tsv',
        rt_perm = expand(
            'analysis/laks_flow/{dataset}/rt_pseudobulks_composite.tsv',
            dataset=[d for d in config['permuted_datasets']]
        )
    output:
        rt_table = 'analysis/laks_flow/permuted/rt_pseudobulks_composite.tsv',
        cor_plot = 'plots/laks_flow/permuted/rt_corr_composite.png',
    params:
        datasets = expand([d for d in config['permuted_datasets']]),

    log: 'logs/laks_flow/permuted/permuted_dataset_rt_profiles.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/permuted_dataset_rt_profiles.py '
        '--rt_ref {input.rt_ref} '
        '--rt_perm {input.rt_perm} '
        '--datasets {params.datasets} '
        '--rt_table {output.rt_table} '
        '--cor_plot {output.cor_plot} '
        '&> {log} ; '
        'deactivate'


rule twidth_analysis_pyro_lf:
    input: 
        cn_T47D = 'analysis/laks_flow/T47D/cn_s_pyro_inferred_filtered.tsv',
        cn_GM18507 = 'analysis/laks_flow/GM18507/cn_s_pyro_inferred_filtered.tsv',
        cn_all = 'analysis/laks_flow/all/cn_s_pyro_inferred_filtered.tsv',
        pseudobulk = 'analysis/laks_flow/all/rt_pseudobulks.tsv'
    output: 
        output_tsv = 'analysis/laks_flow/all/twidth_values.tsv',
        output_png = 'plots/laks_flow/all/twidth_curves.png',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_state = 'model_rep_state'
    log: 'logs/laks_flow/all/twidth_analysis_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_pyro_composite_lf:
    input: 
        cn_T47D = 'analysis/laks_flow/T47D/cn_s_pyro_inferred_composite_filtered.tsv',
        cn_GM18507 = 'analysis/laks_flow/GM18507/cn_s_pyro_inferred_composite_filtered.tsv',
        cn_all = 'analysis/laks_flow/all/cn_s_pyro_inferred_composite_filtered.tsv',
        pseudobulk = 'analysis/laks_flow/all/rt_pseudobulks_composite.tsv'
    output: 
        output_tsv = 'analysis/laks_flow/all/twidth_values_composite.tsv',
        output_png = 'plots/laks_flow/all/twidth_curves_composite.png',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_state = 'model_rep_state'
    log: 'logs/laks_flow/all/twidth_analysis_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_example_cells_lf:
    input:
        cn_t_g1 = 'analysis/laks_flow/T47D/cn_g1.tsv',
        cn_t_g2 = 'analysis/laks_flow/T47D/cn_g2.tsv',
        cn_t_s = 'analysis/laks_flow/T47D/cn_s.tsv',
        cn_gm_g1 = 'analysis/laks_flow/GM18507/cn_g1.tsv',
        cn_gm_g2 = 'analysis/laks_flow/GM18507/cn_g2.tsv',
        cn_gm_s = 'analysis/laks_flow/GM18507/cn_s.tsv',
    output:
        cn_t_g1 = 'plots/laks_flow/T47D/cn_g1_example.png',
        cn_t_g2 = 'plots/laks_flow/T47D/cn_g2_example.png',
        cn_t_s = 'plots/laks_flow/T47D/cn_s_example.png',
        cn_gm_g1 = 'plots/laks_flow/GM18507/cn_g1_example.png',
        cn_gm_g2 = 'plots/laks_flow/GM18507/cn_g2_example.png',
        cn_gm_s = 'plots/laks_flow/GM18507/cn_s_example.png',
    log: 'logs/laks_flow/plot_example_cells.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_example_cells.py '
        '{input} {output} &> {log} ; '
        'deactivate'


rule plot_flow_error_cells_lf:
    input: 'analysis/laks_flow/all/cn_g_pyro_inferred_composite_filtered.tsv'
    output: 'plots/laks_flow/all/flow_error_cells.png'
    log: 'logs/laks_flow/plot_flow_error_cells.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/laks_flow/plot_flow_error_cells.py '
        '{input} {output} &> {log} ; '
        'deactivate'
    