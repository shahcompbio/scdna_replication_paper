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


rule all_sig_lines:
    input:
        expand(
            'plots/sig_lines/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/cn_pseudobulks1.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/twidth_curves.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/clone_rt.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/frac_rt_distributions.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/rpm_embedding.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/subclonal_rt_diffs.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if ((d not in bad_datasets) and (d not in ['OV2295']))
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/phase_changes_confusion.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/signals_heatmaps.pdf',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        'plots/sig_lines/phase_changes_cohort_confusion.png',
        'plots/sig_lines/subclonal_rt_diffs_summary.png',
        'plots/sig_lines/sample_cnas_vs_rt_dists.png',
        'plots/sig_lines/downsampled_twidth_scatter.png',
        'plots/sig_lines/twidth_summary.png',
        'plots/sig_lines/clone_RT_X_profiles.png',
        'plots/sig_lines/clone_corrs.png',
        'plots/sig_lines/frac_rep_distribution.png',
        'plots/sig_lines/Xi_dna_vs_rna_BAF.png'
        
        

def dataset_cn_files_sl(wildcards):
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


def dataset_metric_files_sl(wildcards):
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


def dataset_sample_ids_sl(wildcards):
    if wildcards.dataset == 'OV2295':
        sample_ids = metrics_samples.loc[metrics_samples['isabl_patient_id']==wildcards.dataset]['isabl_sample_id'].unique()
    else:
        col = 'in_{}'.format(wildcards.dataset)
        htert_metrics[col] = list(
            map(lambda x: x.startswith(wildcards.dataset), htert_metrics['isabl_sample_id'])
        )
        sample_ids = htert_metrics.loc[htert_metrics[col]==True]['isabl_sample_id'].unique()
    
    return expand(sample_ids)


rule collect_cn_data_sl:
    input: 
        hmm = dataset_cn_files_sl,
        annotation = dataset_metric_files_sl
    output: 'analysis/sig_lines/{dataset}/cn_data.tsv'
    log: 'logs/sig_lines/{dataset}/collect_cn_data.log'
    params:
        samples = dataset_sample_ids_sl,
        dataset = lambda wildcards: wildcards.dataset
    shell: 
        'python scripts/sig_lines/collect_cn_data.py '
        '--hmm {input.hmm} --annotation {input.annotation} '
        '--samples {params.samples} --dataset {params.dataset} '
        '--output {output} &> {log}'


rule clone_assignments_sl:
    input: 
        cn = 'analysis/sig_lines/{dataset}/cn_data.tsv',
        clones = 'data/signatures/clone_trees/{dataset}_clones.tsv'
    output: 'analysis/sig_lines/{dataset}/cn_data_clones.tsv'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        assign_col = 'copy'
    log: 'logs/sig_lines/{dataset}/clone_assignments.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/clone_assignments.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_ccc_features_sl:
    input: 'analysis/sig_lines/{dataset}/cn_data_clones.tsv'
    output: 'analysis/sig_lines/{dataset}/cn_data_features.tsv'
    log: 'logs/sig_lines/{dataset}/compute_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/compute_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_ccc_features_sl:
    input: 'analysis/sig_lines/{dataset}/cn_data_features.tsv'
    output: 
        plot1 = 'plots/sig_lines/{dataset}/ccc_features_hist.png',
        plot2 = 'plots/sig_lines/{dataset}/ccc_features_scatter.png'
    log: 'logs/sig_lines/{dataset}/plot_ccc_features.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_ccc_features.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule split_cell_cycle_sl:
    input: 'analysis/sig_lines/{dataset}/cn_data_features.tsv'
    output:
        cn_s = 'analysis/sig_lines/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/sig_lines/{dataset}/g1_phase_cells.tsv'
    log: 'logs/sig_lines/{dataset}/split_cell_cycle.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/split_cell_cycle.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule infer_scRT_pyro_sl:
    input:
        cn_s = 'analysis/sig_lines/{dataset}/s_phase_cells.tsv',
        cn_g1 = 'analysis/sig_lines/{dataset}/g1_phase_cells.tsv'
    output:
        main_s_out = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        supp_s_out = 'analysis/sig_lines/{dataset}/scRT_pyro_supp_s_output.tsv',
        main_g_out = 'analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT.tsv',
        supp_g_out = 'analysis/sig_lines/{dataset}/scRT_pyro_supp_g_output.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        cn_prior_method = 'g1_composite',
        infer_mode = 'pyro'
    log: 'logs/sig_lines/{dataset}/infer_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_heatmaps_sl:
    input:
        s_phase = 'analysis/sig_lines/{dataset}/s_phase_cells.tsv',
        g1_phase = 'analysis/sig_lines/{dataset}/g1_phase_cells.tsv'
    output: 'plots/sig_lines/{dataset}/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/sig_lines/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_rt_heatmap_sl:
    input: 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT.tsv'
    output: 'plots/sig_lines/{dataset}/rt_heatmap.png'
    params:
        value_col = 'model_rep_state',
        sort_col = 'model_tau',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_rt_heatmap.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/sig_lines/plot_rt_heatmap.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule plot_pyro_model_output_sl:
    input:
        s_phase = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        g1_phase = 'analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT.tsv'
    output:
        plot1 = 'plots/sig_lines/{dataset}/inferred_cn_rep_results.png',
        plot2 = 'plots/sig_lines/{dataset}/s_vs_g_hmmcopy_states.png',
        plot3 = 'plots/sig_lines/{dataset}/s_vs_g_rpm.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule revise_cell_cycle_labels_sl:
    input: 
        cn_s = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        cn_g = 'analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT.tsv',
    output:
        out_s = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        out_g = 'analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        out_lowqual = 'analysis/sig_lines/{dataset}/model_lowqual_cells.tsv',
    params:
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        rpm_col = 'rpm'
    log: 'logs/sig_lines/{dataset}/revise_cell_cycle_labels.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/revise_cell_cycle_labels.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_filtered_pyro_model_output_sl:
    input:
        s_phase = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g1_phase = 'analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv'
    output:
        plot1 = 'plots/sig_lines/{dataset}/inferred_cn_rep_results_filtered.png',
        plot2 = 'plots/sig_lines/{dataset}/s_vs_g_hmmcopy_states_filtered.png',
        plot3 = 'plots/sig_lines/{dataset}/s_vs_g_rpm_filtered.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_filtered_pyro_model_output.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_pyro_model_output.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT_filtered_sl:
    input: 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 
        plot1 = 'plots/sig_lines/{dataset}/cn_scRT_heatmaps.png',
        plot2 = 'plots/sig_lines/{dataset}/frac_rt_distributions.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'cell_frac_rep'
    log: 'logs/sig_lines/{dataset}/plot_inferred_cn_vs_scRT_filtered.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule cohort_frac_rep_distribution_sl:
    input: 
        expand(
            'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        )
    output: 'plots/sig_lines/frac_rep_distribution.png'
    log: 'logs/sig_lines/cohort_frac_rep_distribution.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/cohort_frac_rep_distribution.py '
        '-i {input} '
        '--plot {output} '
        '&> {log} ; '
        'deactivate'


rule plot_rpm_embedding_sl:
    input:
        s = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        g = 'analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        lowqual = 'analysis/sig_lines/{dataset}/model_lowqual_cells.tsv',
    output: 'plots/sig_lines/{dataset}/rpm_embedding.png'
    params:
        value_col = 'rpm',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_rpm_embedding.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_rpm_embedding.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_rt_pseudobulks_sl:
    input: 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/sig_lines/{dataset}/scRT_pseudobulks.tsv'
    params:
        rep_col = 'model_rep_state',
    log: 'logs/sig_lines/{dataset}/compute_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/compute_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_cn_pseudobulks_sl:
    input: 'analysis/sig_lines/{dataset}/g1_phase_cells.tsv'
    output: 'analysis/sig_lines/{dataset}/cn_pseudobulks.tsv'
    params:
        cn_state_col = 'state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/compute_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/compute_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_cn_pseudobulks_sl:
    input: 'analysis/sig_lines/{dataset}/cn_pseudobulks.tsv'
    output: 
        plot1 = 'plots/sig_lines/{dataset}/cn_pseudobulks1.png',
        plot2 = 'plots/sig_lines/{dataset}/cn_pseudobulks2.png',
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_cn_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_cn_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_clone_rt_pseudobulks_sl:
    input: 
        cn_s = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        rt = 'analysis/sig_lines/{dataset}/scRT_pseudobulks.tsv'
    output: 'plots/sig_lines/{dataset}/clone_rt.png'
    params:
        rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_clone_rt_pseudobulks.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_clone_rt_pseudobulks.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule compute_cell_cycle_clone_counts_sl:
    input:
        cn_s = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        cn_g = 'analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv'
    output: 'analysis/sig_lines/{dataset}/cell_cycle_clone_counts.tsv'
    log: 'logs/sig_lines/{dataset}/compute_cell_cycle_clone_counts.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/compute_cell_cycle_clone_counts.py '
        '{input} {output} &> {log} ; '
        'deactivate'


rule plot_clone_spf_sl:
    input: 'analysis/sig_lines/{dataset}/cell_cycle_clone_counts.tsv'
    output: 'plots/sig_lines/{dataset}/clone_spf.png'
    params:
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/plot_clone_spf.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/plot_clone_spf.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule cohort_clone_counts_sl:
    input:
        expand(
            'analysis/sig_lines/{dataset}/cell_cycle_clone_counts.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
    output: 'analysis/sig_lines/cohort_clone_counts.tsv'
    log: 'logs/sig_lines/cohort_clone_counts.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/cohort_clone_counts.py '
        '--input {input} '
        '--output {output} '
        '&> {log} ; '
        'deactivate'


rule subclonal_rt_diffs_sl:
    input:
        rt = 'analysis/sig_lines/{dataset}/scRT_pseudobulks.tsv',
        cn = 'analysis/sig_lines/{dataset}/cn_pseudobulks.tsv'
    output:
        tsv = 'analysis/sig_lines/{dataset}/subclonal_rt_diffs.tsv',
        png = 'plots/sig_lines/{dataset}/subclonal_rt_diffs.png'
    params:
        rep_col = 'model_rep_state',
        dataset = lambda wildcards: wildcards.dataset
    log: 'logs/sig_lines/{dataset}/subclonal_rt_diffs.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/subclonal_rt_diffs.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule subclonal_rt_diffs_summary_sl:
    input:
        rt = expand(
            'analysis/sig_lines/{dataset}/subclonal_rt_diffs.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in ['OV2295'])
            ]
        )
    output:
        tsv = 'analysis/sig_lines/subclonal_rt_diffs_summary.tsv',
        png = 'plots/sig_lines/subclonal_rt_diffs_summary.png'
    log: 'logs/sig_lines/subclonal_rt_diffs_summary.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/subclonal_rt_diffs_summary.py '
        '-i {input} '
        '--table {output.tsv} '
        '--plot {output.png} '
        '&> {log} ; '
        'deactivate'


rule sample_cnas_vs_rt:
    input:
        rt = expand(
            'analysis/sig_lines/{dataset}/scRT_pseudobulks.tsv',
            dataset=[
                'SA039', 'SA906a', 'SA906b', 'SA1292', 'SA1056', 'SA1188', 'SA1054', 'SA1055'
            ]
        ),
        cn = expand(
            'analysis/sig_lines/{dataset}/cn_pseudobulks.tsv',
            dataset=[
                'SA039', 'SA906a', 'SA906b', 'SA1292', 'SA1056', 'SA1188', 'SA1054', 'SA1055'
            ]
        )
    output:
        plot1 = 'plots/sig_lines/sample_cnas_vs_rt_dists.png',
        plot2 = 'plots/sig_lines/sample_cnas_vs_rt_profiles.png'
    params:
        datasets = expand(['SA039', 'SA906a', 'SA906b', 'SA1292', 'SA1056', 'SA1188', 'SA1054', 'SA1055']),
    log: 'logs/sig_lines/sample_cnas_vs_rt.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/sample_cnas_vs_rt.py '
        '-ir {input.rt} '
        '-ic {input.cn} '
        '-d {params.datasets} '
        '--plot1 {output.plot1} '
        '--plot2 {output.plot2} '
        '&> {log} ; '
        'deactivate'


rule twidth_analysis_sl:
    input: 
        scrt = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        bulks = 'analysis/sig_lines/{dataset}/scRT_pseudobulks.tsv'
    output: 
        output_tsv = 'analysis/sig_lines/{dataset}/twidth_values.tsv',
        output_png = 'plots/sig_lines/{dataset}/twidth_curves.png'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        infer_mode = 'pyro',
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state',
    log: 'logs/sig_lines/{dataset}/twidth_analysis.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/twidth_analysis.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule twidth_summary_sl:
    input: 
        tw = expand(
            'analysis/sig_lines/{dataset}/twidth_values.tsv',
            dataset=[
                'SA039', 'SA906a', 'SA906b', 'SA1292', 'SA1056', 'SA1188', 'SA1054', 'SA1055'
            ]
        )
    output:
        output_tsv = 'analysis/sig_lines/twidth_values.tsv',
        output_png = 'plots/sig_lines/twidth_summary.png'
    params:
        labels = expand(['WT', 'TP53-/-', 'TP53-/-', 'TP53-/-,BRCA1+/-', 'TP53-/-,BRCA1-/-', 'TP53-/-,BRCA2+/-', 'TP53-/-,BRCA2-/-', 'TP53-/-,BRCA2-/-']),
        datasets = expand(['SA039', 'SA906a', 'SA906b', 'SA1292', 'SA1056', 'SA1188', 'SA1054', 'SA1055']),
    log: 'logs/sig_lines/twidth_summary.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/twidth_summary.py '
        '-i {input.tw} '
        '-d {params.datasets} '
        '-l {params.labels} '
        '--table {output.output_tsv} '
        '--plot {output.output_png} '
        '&> {log} ; '
        'deactivate'


rule twidth_downsampling_sl:
    input:
        cn_s = expand(
            'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        )
    output:
        output_tsv = 'analysis/sig_lines/downsampled_twidth_values.tsv',
        output_png = 'plots/sig_lines/downsampled_twidth_scatter.png'
    params:
        datasets = expand([d for d in config['signatures_cell_lines']]),
        labels = expand([str(config['signatures_cell_lines'][d]['type']) for d in config['signatures_cell_lines']]),
        frac_rt_col = 'cell_frac_rep',
        rep_col = 'model_rep_state'
    log: 'logs/sig_lines/twidth_downsampling.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/twidth_downsampling.py '
        '--cn_s {input.cn_s} '
        '-d {params.datasets} '
        '-l {params.labels} '
        '--frac_rt_col {params.frac_rt_col} '
        '--rep_col {params.rep_col} '
        '--table {output.output_tsv} '
        '--plot {output.output_png} '
        '&> {log} ; '
        'deactivate'


rule signals_heatmaps_sl:
    input: 
        ascn = 'analysis/schnapps-results/persample/{dataset}_hscn.csv.gz',
        clones = 'data/signatures/clone_trees/{dataset}_clones.tsv'
    output: 
        figure = 'plots/sig_lines/{dataset}/signals_heatmaps.pdf'
    params:
        dataset = lambda wildcards: wildcards.dataset,
    log: 'logs/sig_lines/{dataset}/signals_heatmaps.log'
    singularity: 'docker://marcjwilliams1/signals'
    # singularity: '/juno/work/shah/users/william1/singularity/signals_v0.7.6.sif'
    shell:
        'Rscript scripts/sig_lines/signals_heatmaps.R '
        '--ascn {input.ascn} '
        '--clones {input.clones} '
        '--dataset {params.dataset} '
        '--heatmap {output.figure} '
        '&> {log}'
    

# TODO: split tsv and png outputs into separate rules
rule phase_changes_sl:
    input:
        cn_g ='analysis/sig_lines/{dataset}/g1_phase_cells_with_scRT_filtered.tsv',
        cn_s ='analysis/sig_lines/{dataset}/s_phase_cells_with_scRT_filtered.tsv',
        cn_lowqual ='analysis/sig_lines/{dataset}/model_lowqual_cells.tsv',
        cn_g_init = 'analysis/sig_lines/{dataset}/g1_phase_cells.tsv',
        cn_s_init = 'analysis/sig_lines/{dataset}/s_phase_cells.tsv',
    output:
        output_tsv = 'analysis/sig_lines/{dataset}/phase_changes.tsv',
        plot1 = 'plots/sig_lines/{dataset}/phase_changes_confusion.png',
        plot2 = 'plots/sig_lines/{dataset}/phase_changes_features.png'
    log: 'logs/sig_lines/{dataset}/phase_changes.log'
    params:
        dataset = lambda wildcards: wildcards.dataset,
        frac_rt_col = 'cell_frac_rep',
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/phase_changes.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule phase_changes_cohort_sl:
    input:
        expand(
            'analysis/sig_lines/{dataset}/phase_changes.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        )
    output:
        plot1 = 'plots/sig_lines/phase_changes_cohort_confusion.png',
        plot2 = 'plots/sig_lines/phase_changes_cohort_features.png'
    log: 'logs/sig_lines/phase_changes_cohort.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/phase_changes_cohort.py '
        '--input {input} '
        '--plot1 {output.plot1} '
        '--plot2 {output.plot2} '
        '&> {log} ; '
        'deactivate'


rule chrX_RT_sl:
    input:
        RT = expand(
            'analysis/sig_lines/{dataset}/scRT_pseudobulks.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        counts = 'analysis/sig_lines/cohort_clone_counts.tsv'
    output:
        sample_rt_profiles = 'plots/sig_lines/sample_RT_X_profiles.png',
        sample_rt_diffs = 'plots/sig_lines/sample_RT_X_diffs.png',
        clone_rt_profiles = 'plots/sig_lines/clone_RT_X_profiles.png',
        clone_rt_diffs_SA1054 = 'plots/sig_lines/SA1054/clone_RT_X_diffs.png',
        clone_rt_diffs_SA1055 = 'plots/sig_lines/SA1055/clone_RT_X_diffs.png',
    log: 'logs/sig_lines/chrX_RT.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/chrX_RT.py '
        '--input {input.RT} '
        '--counts {input.counts} '
        '--sample_rt_profiles {output.sample_rt_profiles} '
        '--sample_rt_diffs {output.sample_rt_diffs} '
        '--clone_rt_profiles {output.clone_rt_profiles} '
        '--clone_rt_diffs_SA1054 {output.clone_rt_diffs_SA1054} '
        '--clone_rt_diffs_SA1055 {output.clone_rt_diffs_SA1055} '
        '&> {log} ; '
        'deactivate'


rule Xi_dna_vs_rna_BAF_sl:
    params:
        datasets = [
            d for d in config['signatures_cell_lines']
            if ((d not in bad_datasets) and (d not in ['SA1292']))
        ]
    output: 
        plot = 'plots/sig_lines/Xi_dna_vs_rna_BAF.png',
        table = 'analysis/sig_lines/Xi_dna_vs_rna_BAF.csv.gz'
    log: 'logs/sig_lines/Xi_dna_vs_rna_BAF.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/Xi_dna_vs_rna_BAF.py '
        '--datasets {params.datasets} '
        '--plot {output.plot} '
        '--table {output.table} '
        '&> {log}'
        ' ; deactivate'


rule cn_and_rt_correlations_sl:
    input:
        rt = expand(
            'analysis/sig_lines/{dataset}/scRT_pseudobulks.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        cn = expand(
            'analysis/sig_lines/{dataset}/cn_pseudobulks.tsv',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets)
            ]
        ),
        counts = 'analysis/sig_lines/cohort_clone_counts.tsv'
    output:
        sample_corrs = 'plots/sig_lines/sample_corrs.png',
        clone_corrs = 'plots/sig_lines/clone_corrs.png',
    log: 'logs/sig_lines/cn_and_rt_correlations.log'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/sig_lines/cn_and_rt_correlations.py '
        '--input_cn {input.cn} '
        '--input_rt {input.rt} '
        '--counts {input.counts} '
        '--sample_corrs {output.sample_corrs} '
        '--clone_corrs {output.clone_corrs} '
        '&> {log} ; '
        'deactivate'
