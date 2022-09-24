import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []
# bad_datasets_3 = ['SA1054', 'SA1055', 'SA1056']
bad_datasets_3 = config['signatures_patient_tumors']

datasets_4 = ['all', 'GM18507', 'T47D']
perm_datasets_4 = [y for x in [datasets_4, config['permuted_datasets']] for y in x]

bad_datasets_5 = config['signatures_cell_lines']

rule all:
    input:
        # fig2 simulated data
        expand(
            'plots/fig2/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/scRT_heatmaps_bulk.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/scRT_heatmaps_pyro_composite.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/twidth_curves_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/twidth_curves_bulk.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/twidth_curves_pyro_composite.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/true_scRT_heatmap.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/cn_vs_scRT_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        expand(
            'plots/fig2/{dataset}/cn_vs_scRT_composite_heatmaps_pyro.png',
            dataset=[
                d for d in config['simulated_datasets']
                if (d not in bad_datasets)
            ]
        ),
        # fig3 hTERT cell lines
        expand(
            'plots/fig3/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_3)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_3)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_3)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_3)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_3)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/inferred_cn_rep_results_nonrep.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_3)
            ]
        ),
        expand(
            'plots/fig3/{dataset}/twidth_curves.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_3)
            ]
        ),
        # fig4 flow-sorted analysis
        expand(
            'plots/fig4/{dataset}/cn_heatmaps.png',
            dataset=[d for d in datasets_4]
        ),
        expand(
            'plots/fig4/{dataset}/cn_clone_heatmaps.png',
            dataset=[d for d in perm_datasets_4]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro.png',
            dataset=[d for d in datasets_4]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro_composite.png',
            dataset=[d for d in perm_datasets_4]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro_filtered.png',
            dataset=[d for d in datasets_4]
        ),
        expand(
            'plots/fig4/{dataset}/scRT_heatmaps_pyro_composite_filtered.png',
            dataset=[d for d in perm_datasets_4]
        ),
        'plots/fig4/all/rt_corr.png',
        'plots/fig4/all/rt_corr_composite.png',
        'plots/fig4/all/twidth_curves.png',
        'plots/fig4/all/twidth_curves_composite.png',
        # fig5 signatures human tumors
        expand(
            'plots/fig5/{dataset}/ccc_features_hist.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/fig5/{dataset}/cn_heatmaps.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/fig5/{dataset}/rt_heatmap.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/fig5/{dataset}/inferred_cn_rep_results.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/fig5/{dataset}/inferred_cn_rep_results_filtered.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/fig5/{dataset}/inferred_cn_rep_results_nonrep.png',
            dataset=[
                d for d in config['signatures_cell_lines']
                if (d not in bad_datasets_5)
            ]
        ),
        expand(
            'plots/fig5/{dataset}/twidth_curves.png',
            dataset=[
                d for d in config['signatures_patient_tumors']
                if (d not in bad_datasets_5)
            ]
        ),

include: "rules/fig2.smk"
include: "rules/fig3.smk"
include: "rules/fig4.smk"
include: "rules/fig5.smk"
