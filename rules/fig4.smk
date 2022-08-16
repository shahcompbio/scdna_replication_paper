import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

training_data = pd.read_csv('/juno/work/shah/users/weinera2/projects/cell_cycle_classifier/cell_cycle_classifier/data/training/curated_feature_data_rt.csv')

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
        'analysis/fig4/cn_data_filtered.tsv',


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

# make sure all cells have same loci
rule filter_data:
    input: 'analysis/fig4/cn_data.tsv'
    output: 'analysis/fig4/cn_data_filtered.tsv'
    log: 'logs/fig4/filter_data.log'
    shell:
        'python scripts/fig4/filter_data.py '
        '{input} {params} {output} &> {log}'
