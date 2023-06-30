import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []

rule all_chrX:
    input:
        'plots/chrX/main_figure.pdf'


rule load_rt_data:
    params:
        rt_paths = config['chrX_rt_paths'],
        count_paths = config['chrX_count_paths']
    output:
        sample_rt = 'analysis/chrX/sample_rt.csv.gz',
        clone_rt = 'analysis/chrX/clone_rt.csv.gz',
        counts = 'analysis/chrX/counts.csv.gz',
    log: 'logs/chrX/load_rt_data.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/chrX/load_rt_data.py '
        '--rt_paths {params.rt_paths} '
        '--count_paths {params.count_paths} '
        '--sample_rt_output {output.sample_rt} '
        '--clone_rt_output {output.clone_rt} '
        '--counts_output {output.counts} '
        '&> {log}'


rule load_signals_data:
    input: 'analysis/chrX/sample_rt.csv.gz',
    output: 'analysis/chrX/sample_arm_bafs.csv.gz',
    log: 'logs/chrX/load_signals_data.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/chrX/load_signals_data.py '
        '--input {input} '
        '--output {output} '
        '&> {log}'
        ' ; deactivate'


rule load_signals_clone_data:
    input: 'analysis/chrX/clone_rt.csv.gz'
    output: 'analysis/chrX/clone_arm_bafs.csv.gz'
    log: 'logs/chrX/load_signals_clone_data.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/chrX/load_signals_clone_data.py '
        '--input {input} '
        '--output {output} '
        '&> {log}'
        ' ; deactivate'


rule plot_chrX_main_figure:
    input:
        sample_bafs = 'analysis/chrX/sample_arm_bafs.csv.gz',
        clone_bafs = 'analysis/chrX/clone_arm_bafs.csv.gz',
        sample_rt = 'analysis/chrX/sample_rt.csv.gz',
        clone_rt = 'analysis/chrX/clone_rt.csv.gz',
        counts = 'analysis/chrX/counts.csv.gz'
    output: 'plots/chrX/main_figure.pdf'
    log: 'logs/chrX/plot_chrX_main_figure.log'
    singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'python3 scripts/chrX/plot_chrX_main_figure.py '
        '--sample_bafs {input.sample_bafs} '
        '--clone_bafs {input.clone_bafs} '
        '--sample_rt {input.sample_rt} '
        '--clone_rt {input.clone_rt} '
        '--counts {input.counts} '
        '--output {output} '
        '&> {log}'