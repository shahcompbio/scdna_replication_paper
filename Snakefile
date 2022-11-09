import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []
# bad_datasets_sl = ['SA1054', 'SA1055', 'SA1056']
bad_datasets_sl = config['signatures_patient_tumors']

datasets_4 = ['all', 'GM18507', 'T47D']
perm_datasets_4 = [y for x in [datasets_4, config['permuted_datasets']] for y in x]

bad_datasets_5 = config['signatures_cell_lines']

rule all:
    input:
        # simulated data
        expand(
            'plots/simulation/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/scRT_heatmaps_bulk.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/scRT_heatmaps_pyro_composite.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/twidth_curves_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/twidth_curves_bulk.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/twidth_curves_pyro_composite.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/true_scRT_heatmap.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/cn_vs_scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/simulation/{dataset}/cn_vs_scRT_composite_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        'plots/simulation/all/model_accuracies.png',
        # sig_lines are hTERT cell lines
        expand(
            'plots/sig_lines/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_sl)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_sl)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_sl)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_sl)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_sl)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/inferred_cn_rep_results_nonrep.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_sl)
            ]
        ),
        expand(
            'plots/sig_lines/{dataset}/twidth_curves.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_sl)
            ]
        ),
        'plots/sig_lines/brca2ko/twidth_curves.png',
        'plots/sig_lines/downsampled_twidth_scatter.png',
        'plots/sig_lines/twidth_summary.png',
        # laks_flow: flow-sorted cell lines from Laks et al
        expand(
            'plots/laks_flow/{dataset}/cn_heatmaps.png',
            dataset=[d for d in datasets_4]
        ),
        expand(
            'plots/laks_flow/{dataset}/cn_clone_heatmaps.png',
            dataset=[d for d in perm_datasets_4]
        ),
        expand(
            'plots/laks_flow/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[d for d in datasets_4]
        ),
        expand(
            'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_composite.png',
            dataset=[d for d in perm_datasets_4]
        ),
        expand(
            'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_filtered.png',
            dataset=[d for d in datasets_4]
        ),
        expand(
            'plots/laks_flow/{dataset}/scRT_heatmaps_pyro_composite_filtered.png',
            dataset=[d for d in perm_datasets_4]
        ),
        expand(
            'analysis/laks_flow/{dataset}/rt_pseudobulks_composite.tsv',
            dataset=[
                d for d in perm_datasets_4
                if (d not in ['T47D', 'GM18507'])
            ]
        ),
        'plots/laks_flow/all/rt_corr.png',
        'plots/laks_flow/all/rt_corr_composite.png',
        'plots/laks_flow/all/twidth_curves.png',
        'plots/laks_flow/all/twidth_curves_composite.png',
        'plots/laks_flow/permuted/summary.png',
        'plots/laks_flow/permuted/rt_corr_composite.png',
        # sig_tumors: signatures human tumors
        expand(
            'plots/sig_tumors/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/inferred_cn_rep_results_nonrep.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/sig_tumors/{dataset}/twidth_curves.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),

include: "rules/simulation.smk"
include: "rules/sig_lines.smk"
include: "rules/laks_flow.smk"
include: "rules/sig_tumors.smk"
