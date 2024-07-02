from argparse import ArgumentParser
import numpy as np
import pandas as pd


def get_args():
    p = ArgumentParser()

    p.add_argument('haplotypes', help='SHAPEIT haplotype phasing of hTERT cell lines')
    p.add_argument('allele_counts', help='table of paths to isabl allele count files')
    p.add_argument('s_phase_cells', help='long-form dataframe of cells in S-phase according to PERT')
    p.add_argument('g_phase_cells', help='long-form dataframe of cells in G1/2-phase according to PERT')
    p.add_argument('dataset', help='name of this dataset')
    p.add_argument('allele_output', help='allele counts but with SIGNALS A/B labels added for each haplotype')
    p.add_argument('s_output', help='table of BAFs for all S-phase cells')
    p.add_argument('g_output', help='table of BAFs for all G1/2-phase cells')

    return p.parse_args()


def load_data(argv, dataset):
    s_phase_cells = pd.read_csv(argv.s_phase_cells, sep='\t', usecols=['cell_id', 'library_id']).drop_duplicates().reset_index(drop=True)
    g_phase_cells = pd.read_csv(argv.g_phase_cells, sep='\t', usecols=['cell_id', 'library_id']).drop_duplicates().reset_index(drop=True)
    allele_counts_table = pd.read_csv(argv.allele_counts)
    
    # compute the set of libraries that are in both the G1-phase and S-phase cells
    g_phase_libraries = set(g_phase_cells['library_id'].unique())
    s_phase_libraries = set(s_phase_cells['library_id'].unique())
    shared_libraries = g_phase_libraries.union(s_phase_libraries)
    print('shared libraries for dataset {}: {}'.format(dataset, shared_libraries))

    temp_table = allele_counts_table.query('dataset == @dataset')
    allele_counts = []
    for i, row in temp_table.iterrows():
        if row['library_id'] in shared_libraries:
            print('loading data for dataset {} library {}'.format(row['dataset'], row['library_id']))
            temp_counts = pd.read_csv(row['allele_counts'], dtype={'chromosome': 'category'})
            print('temp_counts shape: {}'.format(temp_counts.shape))
            allele_counts.append(temp_counts)
    allele_counts = pd.concat(allele_counts, ignore_index=True)

    return s_phase_cells, g_phase_cells, allele_counts


def merge_haplotypes(allele_counts, haplotypes):
    allele_counts['chr'] = allele_counts['chromosome'].astype(str)
    allele_counts['start'] = allele_counts['start'] + 1
    haplotypes['end'] = haplotypes['end'].astype(int)

    # merge the allele counts and haplotypes tables
    merged_df = pd.merge(allele_counts, haplotypes)

    # remove the 'allele' prefix from the phase column and convert to an integer
    merged_df['phase'] = merged_df['phase'].str.replace('allele', '').astype(int)

    # add a new column named 'allele' whose entries are 'A' when the allele_id matches the phase and 'B' otherwise
    merged_df['allele'] = merged_df.apply(lambda x: 'A' if x['allele_id'] == x['phase'] else 'B', axis=1)

    return merged_df


def compute_bafs(s_phase_cells, g_phase_cells, allele_counts):
    s_phase_data = allele_counts[allele_counts['cell_id'].isin(s_phase_cells['cell_id'].values)].groupby(['chr', 'start', 'end', 'hap_label', 'allele'], observed=True)['readcount'].sum().unstack(fill_value=0)
    g_phase_data = allele_counts[allele_counts['cell_id'].isin(g_phase_cells['cell_id'].values)].groupby(['chr', 'start', 'end', 'hap_label', 'allele'], observed=True)['readcount'].sum().unstack(fill_value=0)

    s_phase_data['total'] = s_phase_data[['A', 'B']].sum(axis=1)
    s_phase_data['BAF'] = s_phase_data['B'] / s_phase_data['total']

    g_phase_data['total'] = g_phase_data[['A', 'B']].sum(axis=1)
    g_phase_data['BAF'] = g_phase_data['B'] / g_phase_data['total']

    return s_phase_data, g_phase_data


def main():
    argv = get_args()

    # load in haplotypes
    haplotypes = pd.read_csv(argv.haplotypes)

    # load in allele counts
    s_phase_cells, g_phase_cells, allele_counts = load_data(argv, argv.dataset)

    # merge haplotypes and allele counts
    allele_counts = merge_haplotypes(allele_counts, haplotypes)

    # compute BAFs
    s_phase_data, g_phase_data = compute_bafs(s_phase_cells, g_phase_cells, allele_counts)
    
    # save allele counts
    allele_counts.to_csv(argv.allele_output, index=False)

    # save BAFs
    s_phase_data.to_csv(argv.s_output, index=True)
    g_phase_data.to_csv(argv.g_output, index=True)


if __name__ == '__main__':
    main()
