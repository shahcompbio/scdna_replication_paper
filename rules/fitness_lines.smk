import numpy as np
import pandas as pd

np.random.seed(2794834348)

configfile: "config.yaml"
# samples1 = pd.read_csv('data/fitness/dlp_summaries_rebuttal.csv')
# samples2 = pd.read_csv('data/fitness/fitness_pseudobulk_qc_status.tsv', sep='\t')

# # drop sample_id column from samples2
# samples2.drop(columns=['sample_id'], inplace=True)

# # merge samples1 and samples2 using 'library_id' to get samples
# samples = pd.merge(samples1, samples2, on='library_id')

# # drop rows where ticket, library_id, or sample_id is NA
# samples = samples[samples['mem_jira_ticket'].notna()]
# samples = samples[samples['library_id'].notna()]
# samples = samples[samples['sample_id'].notna()]

bad_datasets = []

rule all_fitness:
    input:
        expand(
            'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv',
            dataset=[
                d for d in config['fitness_lines']
                if (d not in bad_datasets)
            ]
        ),


# use the model output files from sig_lines.smk to assign S-phase cells to timepoints
rule assign_timepoints_fl:
    input: 
        cn_s = 'analysis/sig_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        cn_g = 'analysis/sig_lines/{dataset}/g1_phase_cells.tsv',
        times = 'data/fitness/dlp_summaries_rebuttal.csv'
    output: 
        cn_s = 'analysis/fitness_lines/{dataset}/s_phase_cells_with_scRT.tsv',
        cn_g = 'analysis/fitness_lines/{dataset}/g1_phase_cells.tsv',
    log: 'logs/fitness_lines/{dataset}/assign_timepoints.log'
    shell:
        'python3 scripts/fitness_lines/assign_timepoints.py '
        '{input} {output} &> {log}'
