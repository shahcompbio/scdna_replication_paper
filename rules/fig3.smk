import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"
samples = pd.read_csv('data/signatures/signatures_samples.tsv', sep='\t')

# only look at SA039 and SA906 datasets from fitness paper
# bad_datasets = ['SA1054', 'SA1055', 'SA1056']
bad_datasets = []

rule all_fig3:
    input:
        expand(
            'plots/fig3/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['signatures_datasets']
                if (d not in bad_datasets)
            ]
        ),
        # expand(
        #     'plots/fig3/{dataset}/rt_clusters_heatmap.png',
        #     dataset=[
        #         d for d in config['signatures_datasets']
        #         if (d not in bad_datasets)
        #     ]
        # ),
        # expand(
        #     'plots/fig3/{dataset}/twidth_curves.png',
        #     dataset=[
        #         d for d in config['signatures_datasets']
        #         if (d not in bad_datasets)
        #     ]
        # ),
        expand(
            'plots/fig3/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_datasets']
                if (d not in bad_datasets)
            ]
        ),
        

def dataset_cn_files(wildcards):
    mask = samples['dataset_id'] == wildcards.dataset
    library_ids = samples.loc[mask, 'library_id']
    ticket_ids = samples.loc[mask, 'ticket_id']

    files = expand(
        config['ticket_dir_500kb'] + \
            '/{ticket}/results/hmmcopy/{library}_reads.csv.gz',
        zip, ticket=ticket_ids, library=library_ids
    )
    return files


def dataset_sample_ids(wildcards):
    mask = samples['dataset_id'] == wildcards.dataset
    sample_ids = samples.loc[mask, 'sample_id']
    sample_ids = list(sample_ids)
    return sample_ids


def dataset_metric_files(wildcards):
    mask = samples['dataset_id'] == wildcards.dataset
    library_ids = samples.loc[mask, 'library_id']
    ticket_ids = samples.loc[mask, 'ticket_id']

    files = expand(
        config['ticket_dir_500kb'] + \
            '/{ticket}/results/annotation/{library}_metrics.csv.gz',
        zip, ticket=ticket_ids, library=library_ids
    )
    return files


def dataset_metric_files_updated_classifier(wildcards):
    mask = samples['dataset_id'] == wildcards.dataset
    library_ids = samples.loc[mask, 'library_id']
    sample_ids = samples.loc[mask, 'isabl_sample_id']

    files = expand(
        '/juno/work/shah/users/weinera2/projects/all-dlp-classifier' + \
            '/analysis/{sample}:{library}/merged_classifier_results.tsv',
        zip, sample=sample_ids, library=library_ids
    )
    return files


rule collect_cn_data_3:
    input: 
        hmm = dataset_cn_files,
        annotation = dataset_metric_files_updated_classifier
    output: 'analysis/fig3/{dataset}/cn_data.tsv'
    log: 'logs/fig3/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids
    shell: 
        'python scripts/fig3/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --output {output} &> {log}'


rule clone_assignments_3:
    input: 
        cn = 'analysis/fig3/{dataset}/cn_data.tsv',
        clones ='data/fitness/fitness_cell_assignment_feb07_2020.tsv',
        clones_2295 = 'data/signatures/2295_clones.tsv',
        clones_SA1188 = 'data/signatures/SA1188_clones.tsv',
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
    log: 'logs/fig2/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule get_s_phase_cells_3:
    input: 'analysis/fig3/{dataset}/cn_data_features.tsv'
    output: 'analysis/fig3/{dataset}/s_phase_cells.tsv'
    run:
        df = pd.read_csv(str(input), sep='\t', index_col=False)
        df = df.query('in_tree == False')
        df.to_csv(str(output), sep='\t', index=False)


rule get_non_s_phase_cells_3:
    input: 'analysis/fig3/{dataset}/cn_data_features.tsv'
    output: 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    run:
        df = pd.read_csv(str(input), sep='\t', index_col=False)
        df = df.query('in_tree == True')
        df.to_csv(str(output), sep='\t', index=False)


rule infer_scRT_pyro_3:
    input:
        cn_s = 'analysis/fig3/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_clones',
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
        sort_col = 'model_s_time',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule rt_clustering_3:
    input: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 
        umap = 'plots/fig3/{dataset}/rt_clusters_umap.png',
        kde = 'plots/fig3/{dataset}/rt_clusters_kde.png',
        heatmap = 'plots/fig3/{dataset}/rt_clusters_heatmap.png',
        df = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT_clusters.tsv'
    params:
        value_col = 'rt_state',
        sort_col = 'frac_rt',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/rt_clustering.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/rt_clustering.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule compute_rt_pseudobulks_3:
    input: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'analysis/fig3/{dataset}/scRT_pseudobulks.tsv'
    log: 'logs/fig3/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis_3:
    input: 
        scrt = 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv',
        bulks = 'analysis/fig3/{dataset}/scRT_pseudobulks.tsv'
    output: 'plots/fig3/{dataset}/twidth_curves.png'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/twidth_analysis.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'
