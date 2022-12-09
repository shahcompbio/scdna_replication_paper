import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"
samples1 = pd.read_csv('data/fitness/dlp_summaries_rebuttal.csv')
samples2 = pd.read_csv('data/fitness/fitness_pseudobulk_qc_status.tsv', sep='\t')

# drop sample_id column from samples2
samples2.drop(columns=['sample_id'], inplace=True)

# merge samples1 and samples2 using 'library_id' to get samples
samples = pd.merge(samples1, samples2, on='library_id')

# drop rows where ticket, library_id, or sample_id is NA
samples = samples[samples['mem_jira_ticket'].notna()]
samples = samples[samples['library_id'].notna()]
samples = samples[samples['sample_id'].notna()]

bad_datasets = ['SA609', 'SA000', 'SA039U']

rule all_fitness:
    input:
        expand(
            'plots/fitness/{dataset}/s_vs_g_rpm_filtered.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'analysis/fitness/{dataset}/scRT_pseudobulks.tsv',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/cn_pseudobulks1.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/rpm_embedding.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/ccc_features_scatter.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/clone_spf.png',
            dataset=[
                d for d in config['fitness_rx_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/clone_rt.png',
            dataset=[
                d for d in config['fitness_rx_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fitness/{dataset}/clones_vs_time.png',
            dataset=[
                d for d in config['fitness_datasets']
                if (d not in bad_datasets)
            ]
        ),
        'plots/fitness/fitness_proxy_s_coefficients.png'
        

def dataset_cn_files(wildcards):
    mask = samples['datasetname'] == wildcards.dataset
    library_ids = samples.loc[mask, 'library_id']
    ticket_ids = samples.loc[mask, 'mem_jira_ticket']

    files = expand(
        config['ticket_dir_500kb'] + \
            '/{ticket}/results/hmmcopy/{library}_reads.csv.gz',
        zip, ticket=ticket_ids, library=library_ids
    )
    return files


def dataset_sample_ids(wildcards):
    mask = samples['datasetname'] == wildcards.dataset
    sample_ids = samples.loc[mask, 'sample_id']
    sample_ids = list(sample_ids)
    return sample_ids


def dataset_metric_files(wildcards):
    mask = samples['datasetname'] == wildcards.dataset
    library_ids = samples.loc[mask, 'library_id']
    ticket_ids = samples.loc[mask, 'mem_jira_ticket']

    files = expand(
        config['ticket_dir_500kb'] + \
            '/{ticket}/results/annotation/{library}_metrics.csv.gz',
        zip, ticket=ticket_ids, library=library_ids
    )
    return files


rule collect_cn_data_f:
    input: 
        hmm = dataset_cn_files,
        annotation = dataset_metric_files
    output: 'analysis/fitness/{dataset}/cn_data.tsv'
    log: 'logs/fitness/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids,
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python scripts/common/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --dataset {params.dataset} '
        '--output {output} &> {log}'


rule merge_treatment_groups_f:
    input:
        SA1035U = 'analysis/fitness/SA1035U/cn_data.tsv',
        SA1035T = 'analysis/fitness/SA1035T/cn_data.tsv',
        SA535U = 'analysis/fitness/SA535_CISPLATIN_CombinedU/cn_data.tsv',
        SA535T = 'analysis/fitness/SA535_CISPLATIN_CombinedT/cn_data.tsv'
    output:
        SA1035 = 'analysis/fitness/SA1035/cn_data.tsv',
        SA535 = 'analysis/fitness/SA535/cn_data.tsv'
    log: 'logs/fitness/merge_treatment_groups.log'
    shell:
        'python scripts/fitness/merge_treatment_groups.py '
        '{input} {params} {output} &> {log}'


# make sure all cells have same loci and no NaNs
rule filter_data_f:
    input: 'analysis/fitness/{dataset}/cn_data.tsv'
    output: 'analysis/fitness/{dataset}/cn_data_filtered.tsv'
    log: 'logs/fitness/{dataset}/filter_data.log'
    shell:
        'python scripts/fitness/filter_data.py '
        '{input} {params} {output} &> {log}'


rule assign_timepoints_f:
    input: 
        cn = 'analysis/fitness/{dataset}/cn_data_filtered.tsv',
        times = 'data/fitness/dlp_summaries_rebuttal.csv'
    output: 'analysis/fitness/{dataset}/cn_data_times.tsv'
    log: 'logs/fitness/{dataset}/assign_timepoints.log'
    shell:
        'python3 scripts/fitness/assign_timepoints.py '
        '{input} {output} &> {log}'


rule clone_assignments_f:
    input: 
        cn = 'analysis/fitness/{dataset}/cn_data_times.tsv',
        clones = 'data/fitness/fitness_cell_assignment_feb07_2020.tsv'
    output: 'analysis/fitness/{dataset}/cn_data_clones.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/fitness/{dataset}/clone_assignments.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/clone_assignments.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_ccc_features_f:
    input: 'analysis/fitness/{dataset}/cn_data_clones.tsv'
    output: 'analysis/fitness/{dataset}/cn_data_features.tsv'
    log: 'logs/fitness/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_f:
    input: 'analysis/fitness/{dataset}/cn_data_features.tsv'
    output: 
        plot1 = 'plots/fitness/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/fitness/{dataset}/ccc_features_scatter.png'
    log: 'logs/fitness/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule split_cell_cycle_f:
    input: 'analysis/fitness/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/fitness/{dataset}/s_phase_cells_init_clusters.tsv',
        cn_g1 = 'analysis/fitness/{dataset}/g1_phase_cells_init_clusters.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/split_cell_cycle.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/split_cell_cycle.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule recluster_cells_SA1035_f:
    input: 
        cn_g1 = 'analysis/fitness/SA1035/g1_phase_cells_init_clusters.tsv',
        cn_s = 'analysis/fitness/SA1035/s_phase_cells_init_clusters.tsv'
    output:
        cn_g1 = 'analysis/fitness/SA1035/g1_phase_cells.tsv',
        cn_s = 'analysis/fitness/SA1035/s_phase_cells.tsv',
    params:
        num_clusters = 6
    log: 'logs/fitness/SA1035/recluster_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/recluster_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule recluster_cells_SA535_f:
    input: 
        cn_g1 = 'analysis/fitness/SA535/g1_phase_cells_init_clusters.tsv',
        cn_s = 'analysis/fitness/SA535/s_phase_cells_init_clusters.tsv'
    output:
        cn_g1 = 'analysis/fitness/SA535/g1_phase_cells.tsv',
        cn_s = 'analysis/fitness/SA535/s_phase_cells.tsv',
    params:
        num_clusters = 5
    log: 'logs/fitness/SA535/recluster_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/recluster_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_f:
    input:
        cn_s = 'analysis/fitness/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fitness/{dataset}/g1_phase_cells.tsv'
    output:
        main_s_out = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT.tsv',
        supp_s_out = 'analysis/fitness/{dataset}/scRT_pyro_supp_s_output.tsv',
        main_g_out = 'analysis/fitness/{dataset}/g1_phase_cells_with_scRT.tsv',
        supp_g_out = 'analysis/fitness/{dataset}/scRT_pyro_supp_g_output.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        gc_col = 'gc',
        cn_prior_method = 'g1_clones',
        infer_mode = 'pyro'
    log: 'logs/fitness/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# merge timepoints back in with scRT output
rule merge_scRT_metrics_f:
    input:
        df1 = 'analysis/fitness/{dataset}/s_phase_cells.tsv',
        df2 = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_times.tsv'
    run:
        df1 = pd.read_csv(str(input.df1), sep='\t')  # S-phase cells with timepoints and treatement status mapped to libraries
        df2 = pd.read_csv(str(input.df2), sep='\t')  # scRT model output without timepoints

        # get a mapping of library ID to timepoint & other labels
        lost_metrics = df1[['library_id', 'datasetname', 'label', 'timepoint', 'treated']].drop_duplicates().reset_index(drop=True)

        # merge timepoint back in with scRT model output
        df_out = pd.merge(df2, lost_metrics)

        df_out.to_csv(str(output), sep='\t', index=False)


rule plot_cn_heatmaps_f:
    input:
        s_phase = 'analysis/fitness/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fitness/{dataset}/g1_phase_cells.tsv'
    output: 'plots/fitness/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/common/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_rt_heatmap_f:
    input: 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_times.tsv'
    output: 'plots/fitness/{dataset}/rt_heatmap.png'
    params:
        value_col = 'model_rep_state',
        sort_col = 'model_tau',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/common/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_pyro_model_output_f:
    input:
        s_phase = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_times.tsv',
        g1_phase = 'analysis/fitness/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fitness/{dataset}/inferred_cn_rep_results.png',
        plot2 = 'plots/fitness/{dataset}/s_vs_g_hmmcopy_states.png',
        plot3 = 'plots/fitness/{dataset}/s_vs_g_rpm.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/plot_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/common/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule revise_cell_cycle_labels_f:
    input: 
        cn_s = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_times.tsv',
        cn_g = 'analysis/fitness/{dataset}/g1_phase_cells_with_scRT.tsv',
    output:
        out_s = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        out_g = 'analysis/fitness/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        out_lowqual = 'analysis/fitness/{dataset}/model_lowqual_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/fitness/{dataset}/revise_cell_cycle_labels.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/revise_cell_cycle_labels.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


# split each sample back into treated and untreated groups for plotting and downstream analysis
rule split_by_rx_f:
    input:
        cn_s_SA1035 = 'analysis/fitness/SA1035/s_phase_cells_with_scRT_filtered.tsv',
        cn_g_SA1035 = 'analysis/fitness/SA1035/g1_phase_cells_with_scRT_filtered.tsv',
        cn_s_SA535 = 'analysis/fitness/SA535/s_phase_cells_with_scRT_filtered.tsv',
        cn_g_SA535 = 'analysis/fitness/SA535/g1_phase_cells_with_scRT_filtered.tsv',
    output:
        cn_s_SA1035U = 'analysis/fitness/SA1035U/s_phase_cells_with_scRT_filtered.tsv',
        cn_g_SA1035U = 'analysis/fitness/SA1035U/g1_phase_cells_with_scRT_filtered.tsv',
        cn_s_SA1035T = 'analysis/fitness/SA1035T/s_phase_cells_with_scRT_filtered.tsv',
        cn_g_SA1035T = 'analysis/fitness/SA1035T/g1_phase_cells_with_scRT_filtered.tsv',
        cn_s_SA535U = 'analysis/fitness/SA535_CISPLATIN_CombinedU/s_phase_cells_with_scRT_filtered.tsv',
        cn_g_SA535U = 'analysis/fitness/SA535_CISPLATIN_CombinedU/g1_phase_cells_with_scRT_filtered.tsv',
        cn_s_SA535T = 'analysis/fitness/SA535_CISPLATIN_CombinedT/s_phase_cells_with_scRT_filtered.tsv',
        cn_g_SA535T = 'analysis/fitness/SA535_CISPLATIN_CombinedT/g1_phase_cells_with_scRT_filtered.tsv',
    log: 'logs/fitness/split_by_rx.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/split_by_rx.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_filtered_pyro_model_output_f:
    input:
        s_phase = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/fitness/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
    output:
        plot1 = 'plots/fitness/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/fitness/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/fitness/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/common/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_rpm_embedding_f:
    input:
        s = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g = 'analysis/fitness/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        lowqual = 'analysis/fitness/{dataset}/model_lowqual_cells.tsv',
    output: 'plots/fitness/{dataset}/rpm_embedding.png'
    params:
        value_col = 'rpm',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/plot_rpm_embedding.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/plot_rpm_embedding.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_f:
    input: 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/fitness/{dataset}/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/fitness/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_cn_pseudobulks_f:
    input: 'analysis/fitness/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fitness/{dataset}/cn_pseudobulks.tsv'
    params:
        cn_state_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/compute_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/compute_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_pseudobulks_f:
    input: 'analysis/fitness/{dataset}/cn_pseudobulks.tsv'
    output: 
        plot1 = 'plots/fitness/{dataset}/cn_pseudobulks1.png',
        plot2 = 'plots/fitness/{dataset}/cn_pseudobulks2.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/plot_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/plot_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_clone_rt_and_spf_f:
    input: 
        cn_s = 'analysis/fitness/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        cn_g = 'analysis/fitness/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        rt = 'analysis/fitness/{dataset}/scRT_pseudobulks.tsv'
    output:
        tsv = 'analysis/fitness/{dataset}/cell_cycle_clone_counts.tsv',
        clone_rt = 'plots/fitness/{dataset}/clone_rt.png',
        clone_spf = 'plots/fitness/{dataset}/clone_spf.png'
    params:
        rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fitness/{dataset}/plot_clone_rt_and_spf.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/plot_clone_rt_and_spf.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_clones_vs_time_SA535_f:
    input: 
        treated = 'analysis/fitness/SA535_CISPLATIN_CombinedT/cell_cycle_clone_counts.tsv',
        untreated = 'analysis/fitness/SA535_CISPLATIN_CombinedU/cell_cycle_clone_counts.tsv'
    output: 
        plot1 = 'plots/fitness/SA535/clones_vs_time.png',
        plot2 = 'plots/fitness/SA535/total_cells_vs_time.png'
    params:
        dataset = 'SA535'
    log: 'logs/fitness/SA535/plot_clones_vs_time.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/plot_clones_vs_time.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_clones_vs_time_SA1035_f:
    input: 
        treated = 'analysis/fitness/SA1035T/cell_cycle_clone_counts.tsv',
        untreated = 'analysis/fitness/SA1035U/cell_cycle_clone_counts.tsv'
    output: 
        plot1 = 'plots/fitness/SA1035/clones_vs_time.png',
        plot2 = 'plots/fitness/SA1035/total_cells_vs_time.png'
    params:
        dataset = 'SA1035'
    log: 'logs/fitness/SA1035/plot_clones_vs_time.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/plot_clones_vs_time.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_fitness_proxy_s_coefficients_f:
    input:
        SA1035_treated = 'analysis/fitness/SA1035T/cell_cycle_clone_counts.tsv',
        SA1035_untreated = 'analysis/fitness/SA1035U/cell_cycle_clone_counts.tsv',
        SA535_treated = 'analysis/fitness/SA535_CISPLATIN_CombinedT/cell_cycle_clone_counts.tsv',
        SA535_untreated = 'analysis/fitness/SA535_CISPLATIN_CombinedU/cell_cycle_clone_counts.tsv'
    output: 'plots/fitness/fitness_proxy_s_coefficients.png'
    log: 'logs/fitness/fitness_proxy_s_coefficients.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fitness/plot_fitness_proxy_s_coefficients.py '
        '{input} {output} &> {log} ; '
        'deactivate'


# rule infer_SPF:
#     input:
#         cn_s = 'analysis/fitness/{dataset}/s_phase_cells.tsv',
#         cn_g1 = 'analysis/fitness/{dataset}/non_s_phase_cells.tsv'
#     output:
#         cn_s_out = 'analysis/fitness/{dataset}/s_phase_cells_with_clones.tsv',
#         spf_table = 'analysis/fitness/{dataset}/spf_table.tsv',
#         clone_copy = 'analysis/fitness/{dataset}/clone_copy.tsv'
#     params:
#         input_col = 'copy'
#     log: 'logs/fitness/{dataset}/infer_SPF.log'
#     shell:
#         'source ../scdna_replication_tools/venv/bin/activate ; '
#         'python3 scripts/fitness/infer_SPF.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'


# rule compute_consensus_clone_states:
#     input: 'analysis/fitness/{dataset}/non_s_phase_cells.tsv'
#     output: 'analysis/fitness/{dataset}/clone_states.tsv'
#     params:
#         input_col = 'state',
#         clone_col = 'clone_id'
#     log: 'logs/fitness/{dataset}/compute_consensus_clone_states.log'
#     shell:
#         'source ../scdna_replication_tools/venv/bin/activate ; '
#         'python3 scripts/common/compute_consensus_clone_profiles.py '
#         '{input} {params} {output} &> {log} ; '
#         'deactivate'


# rule plot_consensus_clone_copynumber:
#     input:
#         clone_states = 'analysis/fitness/{dataset}/clone_states.tsv',
#         clone_copy = 'analysis/fitness/{dataset}/clone_copy.tsv'
#     output: 'plots/fitness/{dataset}/consensus_clone_copynumber.pdf'
#     log: 'logs/fitness/{dataset}/plot_consensus_clone_copynumber.log'
#     shell:
#         'source ../scgenome/venv/bin/activate ; '
#         'python3 scripts/fitness/plot_consensus_clone_copynumber.py '
#         '{input} {output} &> {log}'
#         ' ; deactivate'


# rule plot_clone_tree_heatmap:
#     input: 'analysis/fitness/{dataset}/clone_states.tsv',
#     output: 'plots/fitness/{dataset}/clone_tree_heatmap.png'
#     params:
#         dataset = lambda wildcards: wildcards.dataset
#     log: 'logs/fitness/{dataset}/plot_clone_tree_heatmap.log'
#     shell:
#         'source ../scgenome/venv/bin/activate ; '
#         'python3 scripts/fitness/plot_clone_tree_heatmap.py '
#         '{input} {params} {output} &> {log}'
#         ' ; deactivate'


# rule plot_clonal_evolution:
#     input: 
#         s_phase = 'analysis/fitness/{dataset}/s_phase_cells_with_clones.tsv',
#         non_s_phase = 'analysis/fitness/{dataset}/non_s_phase_cells.tsv',
#         times = 'data/fitness/fitness_time_scale.tsv'
#     output:
#         s_out = 'analysis/fitness/{dataset}/s_phase_clone_time_counts.tsv',
#         non_s_out = 'analysis/fitness/{dataset}/non_s_phase_clone_time_counts.tsv',
#         plot = 'plots/fitness/{dataset}/clonal_evolution.pdf'
#     params:
#         dataset = lambda wildcards: wildcards.dataset
#     log: 'logs/fitness/{dataset}/plot_clonal_evolution.log'
#     shell:
#         'python3 scripts/fitness/plot_clonal_evolution.py '
#         '{input} {params} {output} &> {log}'


# rule plot_s_predictiveness:
#     input:
#         s_phase = 'analysis/fitness/{dataset}/s_phase_clone_time_counts.tsv',
#         non_s_phase = 'analysis/fitness/{dataset}/non_s_phase_clone_time_counts.tsv',
#     output:
#         tsv = 'analysis/fitness/{dataset}/s_predictiveness.tsv',
#         plot = 'plots/fitness/{dataset}/s_predictiveness.pdf'
#     params:
#         dataset = lambda wildcards: wildcards.dataset
#     log: 'logs/fitness/{dataset}/plot_s_predictiveness.log'
#     shell:
#         'python3 scripts/fitness/plot_s_predictiveness.py '
#         '{input} {params} {output} &> {log}'


