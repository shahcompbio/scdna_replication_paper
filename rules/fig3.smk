import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"
hmmcopy_samples = pd.read_csv('data/signatures/signatures-hmmcopy.csv')
metrics_samples = pd.read_csv('data/signatures/signatures-annotation.csv')
htert_hmmcopy = hmmcopy_samples.loc[hmmcopy_samples['isabl_patient_id']=='184-hTERT']
htert_metrics = metrics_samples.loc[metrics_samples['isabl_patient_id']=='184-hTERT']

# only look at SA039 and SA906 datasets from fitness paper
# bad_datasets = ['SA1054', 'SA1055', 'SA1056']
bad_datasets = []

rule all_fig3:
    input:
        expand(
            'plots/fig3/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/inferred_cn_rep_results_nonrep.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/twidth_curves.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/clone_rt.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        'plots/fig3/brca2ko/twidth_curves.png',
        'plots/fig3/downsampled_twidth_scatter.png',
        'plots/fig3/twidth_summary.png'
        
        

def dataset_cn_files_3(wildcards):
    if wildcards.dataset == 'OV2295':
        files = hmmcopy_samples.loc[
            hmmcopy_samples['isabl_patient_id']==wildcards.dataset].loc[
            hmmcopy_samples['result_type']=='reads']['result_filepath'].values
    else:
        col = 'in_{}'.format(wildcards.dataset)
        htert_hmmcopy[col] = list(
            map(lambda x: x.startswith(wildcards.dataset), htert_hmmcopy['isabl_sample_id'])
        )
        files = htert_hmmcopy.loc[
            htert_hmmcopy[col]==True].loc[
            htert_hmmcopy['result_type']=='reads']['result_filepath'].values
    
    return expand(files)


def dataset_metric_files_3(wildcards):
    if wildcards.dataset == 'OV2295':
        files = metrics_samples.loc[
            metrics_samples['isabl_patient_id']==wildcards.dataset].loc[
            metrics_samples['result_type']=='metrics']['result_filepath'].values
    else:
        col = 'in_{}'.format(wildcards.dataset)
        htert_metrics[col] = list(
            map(lambda x: x.startswith(wildcards.dataset), htert_metrics['isabl_sample_id'])
        )
        files = htert_metrics.loc[
            htert_metrics[col]==True].loc[
            htert_metrics['result_type']=='metrics']['result_filepath'].values
    
    return expand(files)


def dataset_sample_ids_3(wildcards):
    if wildcards.dataset == 'OV2295':
        sample_ids = metrics_samples.loc[metrics_samples['isabl_patient_id']==wildcards.dataset]['isabl_sample_id'].unique()
    else:
        col = 'in_{}'.format(wildcards.dataset)
        htert_metrics[col] = list(
            map(lambda x: x.startswith(wildcards.dataset), htert_metrics['isabl_sample_id'])
        )
        sample_ids = htert_metrics.loc[htert_metrics[col]==True]['isabl_sample_id'].unique()
    
    return expand(sample_ids)


rule collect_cn_data_3:
    input: 
        hmm = dataset_cn_files_3,
        annotation = dataset_metric_files_3
    output: 'analysis/fig3/{dataset}/cn_data.tsv'
    log: 'logs/fig3/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids_3,
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python scripts/fig3/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --dataset {params.dataset} '
        '--output {output} &> {log}'


rule clone_assignments_3:
    input: 
        cn = 'analysis/fig3/{dataset}/cn_data.tsv',
        clones = 'data/signatures/clone_trees/{dataset}_clones.tsv'
    output: 'analysis/fig3/{dataset}/cn_data_clones.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/fig3/{dataset}/clone_assignments.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/clone_assignments.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_ccc_features_3:
    input: 'analysis/fig3/{dataset}/cn_data_clones.tsv'
    output: 'analysis/fig3/{dataset}/cn_data_features.tsv'
    log: 'logs/fig3/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_3:
    input: 'analysis/fig3/{dataset}/cn_data_features.tsv'
    output: 
        plot1 = 'plots/fig3/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/fig3/{dataset}/ccc_features_scatter.png'
    log: 'logs/fig3/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule split_cell_cycle_3:
    input: 'analysis/fig3/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/fig3/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    log: 'logs/fig3/{dataset}/split_cell_cycle.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/split_cell_cycle.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_3:
    input:
        cn_s = 'analysis/fig3/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output:
        main_out = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv',
        supp_out = 'analysis/fig3/{dataset}/scRT_pyro_supp_output.tsv'  # should contain sample- and library-level params
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pyro'
    log: 'logs/fig3/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_3:
    input:
        s_phase = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv',
        g1_phase = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output: 'plots/fig3/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_rt_heatmap_3:
    input: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'plots/fig3/{dataset}/rt_heatmap.png'
    params:
        value_col = 'model_rep_state',
        sort_col = 'model_tau',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_pyro_model_output_3:
    input:
        s_phase = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv',
        g1_phase = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fig3/{dataset}/inferred_cn_rep_results.png',
        plot2 = 'plots/fig3/{dataset}/s_vs_g_hmmcopy_states.png',
        plot3 = 'plots/fig3/{dataset}/s_vs_g_rpm.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule remove_nonreplicating_cells_3:
    input: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 
        good = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        nonrep = 'analysis/fig3/{dataset}/model_nonrep_cells.tsv',
        lowqual = 'analysis/fig3/{dataset}/model_lowqual_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/fig3/{dataset}/remove_nonreplicating_cells.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/remove_nonreplicating_cells.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_filtered_pyro_model_output_3:
    input:
        s_phase = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fig3/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/fig3/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/fig3/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_nonrep_pyro_model_output_3:
    input:
        s_phase = 'analysis/fig3/{dataset}/model_nonrep_cells.tsv',
        g1_phase = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output:
        plot1 = 'plots/fig3/{dataset}/inferred_cn_rep_results_nonrep.png',
        plot2 = 'plots/fig3/{dataset}/s_vs_g_hmmcopy_states_nonrep.png',
        plot3 = 'plots/fig3/{dataset}/s_vs_g_rpm_nonrep.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_nonrep_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_3:
    input: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/fig3/{dataset}/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/fig3/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_clone_rt_and_spf_3:
    input: 
        cn_s = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        cn_g ='analysis/fig3/{dataset}/g1_phase_cells.tsv',
        cn_g_recovered ='analysis/fig3/{dataset}/model_nonrep_cells.tsv',
        rt = 'analysis/fig3/{dataset}/scRT_pseudobulks.tsv'
    output:
        tsv = 'analysis/fig3/{dataset}/cell_cycle_clone_counts.tsv',
        clone_rt = 'plots/fig3/{dataset}/clone_rt.png',
        clone_spf = 'plots/fig3/{dataset}/clone_spf.png'
    params:
        rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_clone_rt_and_spf.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_clone_rt_and_spf.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_brca2ko_merge_3:
    input: 
        SA1055 = 'analysis/fig3/SA1055/s_phase_cells_with_scRT_filtered.tsv',
        SA1056 = 'analysis/fig3/SA1056/s_phase_cells_with_scRT_filtered.tsv'
    output: 
        cn = 'analysis/fig3/brca2ko/s_phase_cells_with_scRT_filtered.tsv',
        rt_bulks = 'analysis/fig3/brca2ko/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/fig3/brca2ko/compute_rt_pseudobulks_brca2ko_merge.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_rt_pseudobulks_brca2ko_merge.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_3:
    input: 
        scrt = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        bulks = 'analysis/fig3/{dataset}/scRT_pseudobulks.tsv'
    output: 
        output_tsv = 'analysis/fig3/{dataset}/twidth_values.tsv',
        output_png = 'plots/fig3/{dataset}/twidth_curves.png'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        infer_mode = 'pyro',
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
    log: 'logs/fig3/{dataset}/twidth_analysis.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_summary_3:
    input: 
        tw = expand(
            'analysis/fig3/{dataset}/twidth_values.tsv',
            dataset=[
                'SA039', 'SA906a', 'SA906b', 'SA1292', 'SA1056', 'SA1188', 'brca2ko'
            ]
        )
    output:
        output_tsv = 'analysis/fig3/twidth_values.tsv',
        output_png = 'plots/fig3/twidth_summary.png'
    params:
        labels = expand(['WT', 'TP53-/-', 'TP53-/-', 'TP53-/-,BRCA1+/-', 'TP53-/-,BRCA1-/-', 'TP53-/-,BRCA2+/-', 'TP53-/-,BRCA2-/-']),
        datasets = expand(['SA039', 'SA906a', 'SA906b', 'SA1292', 'SA1056', 'SA1188', 'brca2ko']),
    log: 'logs/fig3/twidth_summary.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/twidth_summary.py '
        '-i {input.tw} '
        '-d {params.datasets} '
        '-l {params.labels} '
        '--table {output.output_tsv} '
        '--plot {output.output_png} '
        '&> {log} ; '
        'deactivate'


rule twidth_downsampling_3:
    input:
        cn_s = expand(
            'analysis/fig3/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        )
    output:
        output_tsv = 'analysis/fig3/downsampled_twidth_values.tsv',
        output_png = 'plots/fig3/downsampled_twidth_scatter.png'
    params:
        datasets = expand([d for d in config['signatures_cell_lines']]),
        labels = expand([str(config['signatures_cell_lines'][d]['type']) for d in config['signatures_cell_lines']]),
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state'
    log: 'logs/fig3/twidth_downsampling.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/twidth_downsampling.py '
        '--cn_s {input.cn_s} '
        '-d {params.datasets} '
        '-l {params.labels} '
        '--frac_rt_col {params.frac_rt_col} '
        '--rep_col {params.rep_col} '
        '--table {output.output_tsv} '
        '--plot {output.output_png} '
        '&> {log} ; '
        'deactivate'