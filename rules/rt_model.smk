import pandas as pd
import numpy as np

np.random.seed(2794834348)

configfile: "config.yaml"

bad_datasets = []

rule all_rt_model:
    input:
        'plots/rt_model/loss_curve.pdf',
        'plots/rt_model/loss_curve_noX.pdf',
        'plots/rt_model/beta_importance_posteriors.pdf',
        'plots/rt_model/beta_importance_posteriors_noX.pdf',
        'plots/rt_model/rt_profile_posteriors.pdf',
        'plots/rt_model/rt_profile_posteriors_noX.pdf'


rule load_data_rt:
    params:
        table = 'data/rt_model_metacohort.tsv',
        samples = config['rt_model_datasets'],
        count_paths = config['rt_model_count_paths']
    output:
        sample_rt = 'analysis/rt_model/sample_rt.csv.gz',
        clone_rt = 'analysis/rt_model/clone_rt.csv.gz',
        sample_cn = 'analysis/rt_model/sample_cn.csv.gz',
        clone_cn = 'analysis/rt_model/clone_cn.csv.gz',
        features = 'analysis/rt_model/features.csv.gz',
    log: 'logs/rt_model/load_data.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/load_data.py '
        '--table {params.table} '
        '--samples {params.samples} '
        '--count_paths {params.count_paths} '
        '--sample_rt {output.sample_rt} '
        '--clone_rt {output.clone_rt} '
        '--sample_cn {output.sample_cn} '
        '--clone_cn {output.clone_cn} '
        '--features {output.features} '
        '&> {log}'
        ' ; deactivate'


rule run_rt_model:
    input:
        clone_rt = 'analysis/rt_model/clone_rt.csv.gz',
        features = 'analysis/rt_model/features.csv.gz'
    output:
        beta_importance_posteriors = 'analysis/rt_model/beta_importance_posteriors.csv.gz',
        rt_profile_posteriors = 'analysis/rt_model/rt_profile_posteriors.csv.gz',
        loss_curve = 'analysis/rt_model/loss_curve.csv.gz',
    log: 'logs/rt_model/run_rt_model.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/rt_model.py '
        '--clone_rt {input.clone_rt} '
        '--features {input.features} '
        '--beta_importance_posteriors {output.beta_importance_posteriors} '
        '--rt_profile_posteriors {output.rt_profile_posteriors} '
        '--loss_curve {output.loss_curve} '
        '&> {log}'
        ' ; deactivate'


rule run_rt_model_noX:
    input:
        clone_rt = 'analysis/rt_model/clone_rt.csv.gz',
        features = 'analysis/rt_model/features.csv.gz'
    output:
        beta_importance_posteriors = 'analysis/rt_model/beta_importance_posteriors_noX.csv.gz',
        rt_profile_posteriors = 'analysis/rt_model/rt_profile_posteriors_noX.csv.gz',
        loss_curve = 'analysis/rt_model/loss_curve_noX.csv.gz',
    log: 'logs/rt_model/run_rt_model_noX.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/rt_model.py '
        '--clone_rt {input.clone_rt} '
        '--features {input.features} '
        '--remove_x '
        '--beta_importance_posteriors {output.beta_importance_posteriors} '
        '--rt_profile_posteriors {output.rt_profile_posteriors} '
        '--loss_curve {output.loss_curve} '
        '&> {log}'
        ' ; deactivate'


rule plot_loss_curve_rt:
    input: 'analysis/rt_model/loss_curve.csv.gz'
    output: 'plots/rt_model/loss_curve.pdf'
    log: 'logs/rt_model/plot_loss_curve.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/plot_loss_curve.py '
        '{input} {output} &> {log}'
        ' ; deactivate'


rule plot_loss_curve_noX_rt:
    input: 'analysis/rt_model/loss_curve_noX.csv.gz'
    output: 'plots/rt_model/loss_curve_noX.pdf'
    log: 'logs/rt_model/plot_loss_curve_noX.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/plot_loss_curve.py '
        '{input} {output} &> {log}'
        ' ; deactivate'


rule plot_beta_importance_posteriors_rt:
    input: 'analysis/rt_model/beta_importance_posteriors.csv.gz'
    output: 'plots/rt_model/beta_importance_posteriors.pdf'
    log: 'logs/rt_model/plot_beta_importance_posteriors.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/plot_beta_importance_posteriors.py '
        '{input} {output} &> {log}'
        ' ; deactivate'


rule plot_beta_importance_posteriors_noX_rt:
    input: 'analysis/rt_model/beta_importance_posteriors_noX.csv.gz'
    output: 'plots/rt_model/beta_importance_posteriors_noX.pdf'
    log: 'logs/rt_model/plot_beta_importance_posteriors_noX.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/plot_beta_importance_posteriors.py '
        '{input} --noX {output} &> {log}'
        ' ; deactivate'


rule plot_rt_profile_posteriors_rt:
    input: 'analysis/rt_model/rt_profile_posteriors.csv.gz'
    output: 'plots/rt_model/rt_profile_posteriors.pdf'
    log: 'logs/rt_model/plot_rt_profile_posteriors.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/plot_rt_profile_posteriors.py '
        '{input} {output} &> {log}'
        ' ; deactivate'


rule plot_rt_profile_posteriors_noX_rt:
    input: 'analysis/rt_model/rt_profile_posteriors_noX.csv.gz'
    output: 'plots/rt_model/rt_profile_posteriors_noX.pdf'
    log: 'logs/rt_model/plot_rt_profile_posteriors_noX.log'
    # singularity: 'docker://adamcweiner/scdna_replication_tools:main'
    shell:
        'source ../scdna_replication_tools/venv3/bin/activate ; '
        'python3 scripts/rt_model/plot_rt_profile_posteriors.py '
        '{input} --noX {output} &> {log}'
        ' ; deactivate'