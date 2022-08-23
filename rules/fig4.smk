import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

training_data = pd.read_csv('/juno/work/shah/users/weinera2/projects/cell_cycle_classifier/cell_cycle_classifier/data/training/curated_feature_data_rt.csv')

# paths on juno for reads, metics, and annotation metrics
cn_data_urls = [
    '/work/shah/tantalus/SC-1563/results/results/hmmcopy_autoploidy/A90553C_multiplier0_reads.csv.gz',
    '/work/shah/tantalus/SC-1561/results/results/hmmcopy_autoploidy/A73044A_multiplier0_reads.csv.gz',
    '/work/shah/tantalus/SC-1583/results/results/hmmcopy_autoploidy/A96139A_multiplier0_reads.csv.gz',
    '/work/shah/tantalus/SC-1585/results/results/hmmcopy_autoploidy/A96147A_multiplier0_reads.csv.gz',
]

metrics_data_urls = [
    '/work/shah/tantalus/SC-1563/results/results/hmmcopy_autoploidy/A90553C_multiplier0_metrics.csv.gz',
    '/work/shah/tantalus/SC-1561/results/results/hmmcopy_autoploidy/A73044A_multiplier0_metrics.csv.gz',
    '/work/shah/tantalus/SC-1583/results/results/hmmcopy_autoploidy/A96139A_multiplier0_metrics.csv.gz',
    '/work/shah/tantalus/SC-1585/results/results/hmmcopy_autoploidy/A96147A_multiplier0_metrics.csv.gz',
]

align_metrics_data_urls = [
    '/work/shah/tantalus/SC-1563/results/results/alignment/A90553C_alignment_metrics.csv.gz',
    '/work/shah/tantalus/SC-1561/results/results/alignment/A73044A_alignment_metrics.csv.gz',
    '/work/shah/tantalus/SC-1583/results/results/alignment/A96139A_alignment_metrics.csv.gz',
    '/work/shah/tantalus/SC-1585/results/results/alignment/A96147A_alignment_metrics.csv.gz',
]

rule all_fig4:
    input:
        'analysis/fig4/cn_s_pyro_infered.tsv',
        'plots/fig4/cn_heatmaps.png',
        'plots/fig4/scRT_heatmaps_pyro.png',


# fetch the raw data
rule get_data:
    input:
        cn_data_urls = cn_data_urls,
        metrics_data_urls = metrics_data_urls,
        align_metrics_data_urls = align_metrics_data_urls,
    output:
        cn_data = 'analysis/fig4/cn_data.tsv',
        metrics_data = 'analysis/fig4/metrics_data.tsv',
        align_metrics_data = 'analysis/fig4/align_metrics_data.tsv'
    log: 'logs/fig4/get_data.log'
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
rule filter_data:
    input: 'analysis/fig4/cn_data.tsv'
    output: 'analysis/fig4/cn_data_filtered.tsv'
    log: 'logs/fig4/filter_data.log'
    shell:
        'python scripts/fig4/filter_data.py '
        '{input} {params} {output} &> {log}'


# use metrics file to split each cell in filtered cn data by cell cycle state
rule split_cell_cycle:
    input:
        cn_data = 'analysis/fig4/cn_data_filtered.tsv',
        metrics_data = 'analysis/fig4/metrics_data.tsv'
    output:
        cn_s = 'analysis/fig4/cn_s.tsv',
        cn_g1 = 'analysis/fig4/cn_g1.tsv',
        cn_g2 = 'analysis/fig4/cn_g2.tsv'
    log: 'logs/fig4/split_cell_cycle.log'
    shell:
        'python scripts/fig4/split_cell_cycle.py '
        '{input} {params} {output} &> {log}'


rule plot_cn_heatmaps:
    input:
        s_phase = 'analysis/fig4/cn_s.tsv',
        g1_phase = 'analysis/fig4/cn_g1.tsv',
        g2_phase = 'analysis/fig4/cn_g2.tsv',
    output: 'plots/fig4/cn_heatmaps.png'
    params:
        value_col = 'state',
        dataset = 'flow-sorted'
    log:
        'logs/fig4/plot_cn_heatmaps.log'
    shell:
        'source ../scgenome/venv/bin/activate ; '
        'python3 scripts/fig4/plot_s_vs_g_cn_heatmaps.py '
        '{input} {params} {output} &> {log}'
        ' ; deactivate'


rule compute_cn_prior:
    input:
        cn_s = 'analysis/fig4/cn_s.tsv',
        cn_g1 = 'analysis/fig4/cn_g1.tsv',
        cn_g2 = 'analysis/fig4/cn_g2.tsv'
    output: 'analysis/fig4/cn_s_with_cn_prior.tsv',
    log: 'logs/fig4/compute_cn_prior.log'
    shell:
        'python scripts/fig4/compute_cn_prior.py '
        '{input} {params} {output} &> {log}'


rule infer_scRT_pyro:
    input:
        cn_s = 'analysis/fig4/cn_s.tsv',
        cn_g1 = 'analysis/fig4/cn_g1.tsv',
        cn_g2 = 'analysis/fig4/cn_g2.tsv'
    output: 'analysis/fig4/cn_s_pyro_infered.tsv',
    params:
        input_col = 'rpm',
        cn_col = 'state',
        copy_col = 'copy',
        gc_col = 'gc',
        infer_mode = 'pyro'
    log: 'logs/fig4/infer_scRT_pyro.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/infer_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'


rule plot_inferred_cn_vs_scRT:
    input: 'analysis/fig4/cn_s_pyro_infered.tsv'
    output: 
        plot1 = 'plots/fig4/scRT_heatmaps_pyro.png',
        plot2 = 'plots/fig4/frac_rt_distributions_pyro.png'
    params:
        rep_col = 'model_rep_state',
        cn_col = 'model_cn_state',
        frac_rt_col = 'model_s_time'
    log: 'logs/fig4/plot_inferred_cn_vs_scRT.log'
    shell:
        'source ../scdna_replication_tools/venv/bin/activate ; '
        'python3 scripts/fig4/plot_inferred_cn_vs_scRT.py '
        '{input} {params} {output} &> {log} ; '
        'deactivate'
