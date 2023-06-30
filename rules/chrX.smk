import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []

rule all_chrX:
    input:
        'analysis/chrX/sample_arm_bafs.csv.gz',


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
