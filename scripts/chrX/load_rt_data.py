import pandas as pd
import numpy as np
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-rp', '--rt_paths', type=str, nargs='+', help='list of pseudobulk RT profiles for each sample')
    p.add_argument('-cp', '--count_paths', type=str, nargs='+', help='cell cycle clone counts for across all samples (split across multiple files)')
    p.add_argument('-sro', '--sample_rt_output', type=str, help='chrX RT delays for sample profiles')
    p.add_argument('-cro', '--clone_rt_output', type=str, help='chrX RT delays for clone profiles')
    p.add_argument('-co', '--counts_output', type=str, help='merged counts for all samples (merged into one file)')

    return p.parse_args()


def load_rt(argv):
    # read in the pseudobulk RT profiles for each sample
    rt = pd.DataFrame()
    for path in argv.rt_paths:
        if path.endswith('tsv'):
            temp_rt = pd.read_csv(path, sep='\t')
        elif path.endswith('csv.gz'):
            temp_rt = pd.read_csv(path)
        cols = [c for c in temp_rt.columns if c.startswith('pseudobulk') and c.endswith('model_rep_state')]
        cols.append('chr')
        cols.append('start')
        temp_rt = temp_rt[cols]
        
        # add dataset name as prefix to RT columns
        d = path.split('/')[-2]
        for c in temp_rt.columns:
            if c.startswith('pseudobulk') and c.endswith('model_rep_state'):
                temp_rt.rename(columns={c: '{}_{}'.format(d, c)}, inplace=True)
        
        if rt.empty:
            rt = temp_rt
        else:
            rt = pd.merge(rt, temp_rt)

    # set chr column to category
    rt.chr = rt.chr.astype('str')
    rt.chr = rt.chr.astype('category')

    # add end position as it's necessary for plotting functions
    rt['end'] = rt['start'] + 500000 - 1

    return rt


def load_counts(argv):
    counts = []
    # load counts for each file provided
    for path in argv.count_paths:
        if path.endswith('tsv'):
            temp_counts = pd.read_csv(path, sep='\t')
        elif path.endswith('csv.gz'):
            temp_counts = pd.read_csv(path)
        # append this file's counts to the list
        counts.append(temp_counts)
    # concatenate all counts into one dataframe
    counts = pd.concat(counts, ignore_index=True)

    return counts


def compute_mean_chrX_rt_delay(rt):
    ''' For each sample or clone RT profile, find the mean RT difference between autosomes and chrX. '''
    # loop through all columns in the RT dataframe that aren't chr, start, or end
    mean_rt = []
    for c in rt.columns:
        if c not in ['chr', 'start', 'end']:
            # find that average RT value for autosomes (chr1-22)
            autosomes = np.mean(rt.query('chr!="X" and chr!="Y"')[c].values)
            autosomes_std = np.std(rt.query('chr!="X" and chr!="Y"')[c].values)
            # find the average RT value for chrX
            chrXp = np.mean(rt.query('chr=="X"').query('start<60000000')[c].values)
            chrXq = np.mean(rt.query('chr=="X"').query('start>60000000')[c].values)
            chrX = np.mean(rt.query('chr=="X"')[c].values)
            # find the standard deviation RT for chrX
            chrXp_std = np.std(rt.query('chr=="X"').query('start<60000000')[c].values)
            chrXq_std = np.std(rt.query('chr=="X"').query('start>60000000')[c].values)
            chrX_std = np.std(rt.query('chr=="X"')[c].values)
            # find the dataset that this column belongs to
            dataset = c.split('_')[0]
            # find out whether this column is clone or sample RT profile
            clone_id = 'sample'
            if 'clone' in c:
                clone_id = c.split('_')[2].replace('clone', '')
            # create a temporary dataframe for this column
            temp_df = pd.DataFrame({'dataset': [dataset], 'clone_id': [clone_id], 
                'autosome_mean_rt': [autosomes], 'autosome_std_rt': [autosomes_std],
                'chrXp_mean_rt': [chrXp], 'chrXq_mean_rt': [chrXq], 'chrX_mean_rt': [chrX],
                'chrXp_std_rt': [chrXp_std], 'chrXq_std_rt': [chrXq_std], 'chrX_std_rt': [chrX_std]
                })
            mean_rt.append(temp_df)
    mean_rt = pd.concat(mean_rt, ignore_index=True)

    # create a new column named 'mean_chrX_rt_delay' that is the difference between autosome and chrX mean RT values
    mean_rt['mean_chrXp_rt_delay'] = mean_rt['chrXp_mean_rt'] - mean_rt['autosome_mean_rt']
    mean_rt['mean_chrXq_rt_delay'] = mean_rt['chrXq_mean_rt'] - mean_rt['autosome_mean_rt']
    mean_rt['mean_chrX_rt_delay'] = mean_rt['chrX_mean_rt'] - mean_rt['autosome_mean_rt']
    # standard deviations need to be added together to get uncertainty estimates for the mean RT delay
    mean_rt['std_chrXp_rt_delay'] = mean_rt['chrXp_std_rt'] + mean_rt['autosome_std_rt']
    mean_rt['std_chrXq_rt_delay'] = mean_rt['chrXq_std_rt'] + mean_rt['autosome_std_rt']
    mean_rt['std_chrX_rt_delay'] = mean_rt['chrX_std_rt'] + mean_rt['autosome_std_rt']

    return mean_rt


def main():
    argv = get_args()

    # load RT profiles
    rt = load_rt(argv)

    # load counts
    counts = load_counts(argv)

    # compute mean RT delay for each sample and clone
    mean_rt = compute_mean_chrX_rt_delay(rt)

    # split mean_rt into sample and clone dataframes
    sample_mean_rt = mean_rt.query('clone_id=="sample"')
    clone_mean_rt = mean_rt.query('clone_id!="sample"')

    # save the output files
    sample_mean_rt.to_csv(argv.sample_rt_output, index=False)
    clone_mean_rt.to_csv(argv.clone_rt_output, index=False)
    counts.to_csv(argv.counts_output, index=False)


if __name__ == '__main__':
    main()
