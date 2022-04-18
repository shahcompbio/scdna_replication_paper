import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"
samples = pd.read_csv('data/signatures/signatures_samples.tsv', sep='\t')

# only look at SA039 and SA906 datasets from fitness paper
bad_datasets = ['SA1054', 'SA1055', 'SA1056']

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
        expand(
            'plots/fig3/{dataset}/rt_clusters_heatmap.png',
            dataset=[
                d for d in config['signatures_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/twidth_curves.png',
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


rule collect_cn_data:
    input: 
        hmm = dataset_cn_files,
        annotation = dataset_metric_files_updated_classifier
    output: 'analysis/fig3/{dataset}/{dataset}_cn_data.tsv'
    log: 'logs/fig3/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids
    shell: 
        'python scripts/fig3/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --output {output} &> {log}'


rule get_s_phase_cells:
    input: 'analysis/fig3/{dataset}/{dataset}_cn_data.tsv'
    output: 'analysis/fig3/{dataset}/s_phase_cells.tsv'
    run:
        df = pd.read_csv(str(input), sep='\t', index_col=False)
        df = df.query('is_s_phase_prob_new >= 0.5 | is_s_phase_prob >= 0.5')
        df.to_csv(str(output), sep='\t', index=False)


rule get_non_s_phase_cells:
    input:
        cn = 'analysis/fig3/{dataset}/{dataset}_cn_data.tsv',
        clones = 'data/fitness/fitness_cell_assignment_feb07_2020.tsv'
    output: 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset
    run:
        df = pd.read_csv(str(input.cn), sep='\t', index_col=False)
        #clones = pd.read_csv(str(input.clones), sep='\t', index_col=False)
        dataset = str(params.dataset)

        if dataset in ['2295', 'SA1188']:
            clones = pd.read_csv('data/signatures/{}_clones.tsv'.format(dataset), sep='\t')
        else:
            # load in clones from fitness results
            clones = pd.read_csv(str(input.clones))
            clones = clones.drop(columns=['V1', 'datatag', 'sample_id'])
            clones = clones.rename(columns={'single_cell_id': 'cell_id', 'letters': 'clone_id'})

            # remove the 'a' or 'b' suffix from SA906 cell IDs in the clone mapping file
            if 'SA906' in dataset:
                clones['cell_id'] = clones['cell_id'].apply(lambda x: x.replace(dataset, 'SA906'))

        df = df.query('is_s_phase_prob_new < 0.5 & is_s_phase_prob < 0.5')
        # only use cells that have clone_id's assigned (and add clone_id column)
        df = pd.merge(df, clones, on='cell_id')

        df.to_csv(str(output), sep='\t', index=False)


rule infer_scRT:
    input:
        cn_s = 'analysis/fig3/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv',
    params:
        input_col = 'reads',
        assign_col = 'copy',
        infer_mode = 'cell'
    log: 'logs/fig3/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_g1:
    input:
        cn_s = 'analysis/fig3/{dataset}/g1_phase_cells.tsv',
        cn_g1 = 'analysis/fig3/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/fig3/{dataset}/g1_phase_cells_with_scRT.tsv',
    params:
        input_col = 'reads',
        assign_col = 'copy',
        infer_mode = 'clone'
    log: 'logs/fig3/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps:
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


rule plot_rt_heatmap:
    input: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'plots/fig3/{dataset}/rt_heatmap.png'
    params:
        value_col = 'rt_state',
        sort_col = 'frac_rt',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/fig3/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig3/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule rt_clustering:
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


rule compute_rt_pseudobulks:
    input: 'analysis/fig3/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'analysis/fig3/{dataset}/scRT_pseudobulks.tsv'
    log: 'logs/fig3/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig3/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_analysis:
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
