import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"
hmmcopy_samples = pd.read_csv('data/signatures/signatures-hmmcopy.csv')
metrics_samples = pd.read_csv('data/signatures/signatures-annotation.csv')

bad_datasets = []

rule all_fig5:
    input:
        expand(
            'plots/fig5/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        # expand(
        #     'plots/fig5/{dataset}/rt_heatmap.png',
        #     dataset=[
        #         d for d in config['signatures_patient_tumors']
        #         if (d not in bad_datasets)
        #     ]
        # ),
        expand(
            'plots/fig5/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets)
            ]
        ),
        

def dataset_cn_files_5(wildcards):
    files = hmmcopy_samples.loc[
        hmmcopy_samples['isabl_patient_id']==wildcards.dataset].loc[
        hmmcopy_samples['result_type']=='reads']['result_filepath'].values
    return expand(files)


def dataset_metric_files_5(wildcards):
    files = metrics_samples.loc[
        metrics_samples['isabl_patient_id']==wildcards.dataset].loc[
        metrics_samples['result_type']=='metrics']['result_filepath'].values
    return expand(files)


def dataset_sample_ids_5(wildcards):
    sample_ids = metrics_samples.loc[metrics_samples['isabl_patient_id']==wildcards.dataset]['isabl_sample_id'].unique()
    return expand(sample_ids)


rule collect_cn_data_5:
    input: 
        hmm = dataset_cn_files_5,
        annotation = dataset_metric_files_5
    output: 'analysis/fig5/{dataset}/cn_data.tsv'
    log: 'logs/fig5/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids_5,
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python scripts/fig5/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --dataset {params.dataset} '
        '--output {output} &> {log}'


rule clone_assignments_5:
    input: 
        cn = 'analysis/fig5/{dataset}/cn_data.tsv',
        clones = 'data/signatures/clone_trees/{dataset}_clones.tsv'
    output: 'analysis/fig5/{dataset}/cn_data_clones.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/fig5/{dataset}/clone_assignments.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/clone_assignments.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_ccc_features_5:
    input: 'analysis/fig5/{dataset}/cn_data_clones.tsv'
    output: 'analysis/fig5/{dataset}/cn_data_features.tsv'
    log: 'logs/fig5/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_5:
    input: 'analysis/fig5/{dataset}/cn_data_features.tsv'
    output: 
        plot1 = 'plots/fig5/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/fig5/{dataset}/ccc_features_scatter.png'
    log: 'logs/fig5/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig5/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule get_s_phase_cells_5:
    input: 'analysis/fig5/{dataset}/cn_data_features.tsv'
    output: 'analysis/fig5/{dataset}/s_phase_cells.tsv'
    run:
        df = pd.read_csv(str(input), sep='\t', index_col=False)
        df = df.query('in_tree == False')
        df.to_csv(str(output), sep='\t', index=False)


rule get_non_s_phase_cells_5:
    input: 'analysis/fig5/{dataset}/cn_data_features.tsv'
    output: 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    run:
        df = pd.read_csv(str(input), sep='\t', index_col=False)
        df = df.query('in_tree == True')
        df.to_csv(str(output), sep='\t', index=False)


rule infer_scRT_pyro_5:
    input:
        cn_s = 'analysis/fig5/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fig5/{dataset}/s_phase_cells_with_scRT.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_clones',
        infer_mode = 'pyro'
    log: 'logs/fig5/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_5:
    input:
        s_phase = 'analysis/fig5/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/fig5/{dataset}/g1_phase_cells.tsv'
    output: 'plots/fig5/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig5/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_rt_heatmap_5:
    input: 'analysis/fig5/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'plots/fig5/{dataset}/rt_heatmap.png'
    params:
        value_col = 'model_rep_state',
        sort_col = 'model_s_time',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig5/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'
