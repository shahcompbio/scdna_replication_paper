from argparse import ArgumentParser
import pandas as pd


def get_args():
    p = ArgumentParser()

    p.add_argument('--hmm', type=str, nargs='+', help='list of hmm reads data from different library and ticket ids')
    p.add_argument('--annotation', type=str, nargs='+', help='list of annotation metrics data from different library and ticket ids')
    p.add_argument('--samples', type=str, nargs='+', help='list of samples that should be included, all others (normal or controls) will be filtered out')
    p.add_argument('--dataset', type=str, help='name of this dataset -- used for filtering out control cells')
    p.add_argument('--output', type=str, help='output tsv file that combines all cells from same dataset into one file')

    return p.parse_args()


def compute_reads_per_million(cn, reads_col='reads', rpm_col='rpm'):
    for cell_id, cell_cn in cn.groupby('cell_id'):
        x = cell_cn[reads_col].values
        temp_rpm = x * 1e6 / sum(x)
        cn.loc[cell_cn.index, rpm_col] = temp_rpm
    return cn


if __name__ == '__main__':
    argv = get_args()
    cn_pieces = []

    # load the hmmcopy reads data
    for f in argv.hmm:
        piece = pd.read_csv(
            f, index_col=['chr', 'start', 'end', 'cell_id'],
            # dtype=str
        )
        piece = piece[['reads', 'gc', 'map', 'copy', 'state']]
        cn_pieces.append(piece)

    cn = pd.concat(cn_pieces)
    cn = cn.reset_index()

    # load the hmmcopy annotation data
    met_pieces = []
    for f in argv.annotation:
        piece = pd.read_csv(
            f, index_col=['cell_id'],
        )
        piece = piece[[
            'total_mapped_reads_hmmcopy', 'experimental_condition', 'is_s_phase', 'is_s_phase_prob', 'multiplier',
            'quality', 'is_contaminated', 'cell_call', 'fastqscreen_grch37', 'fastqscreen_mm10', 'coverage_depth',
            'breakpoints'
        ]]
        met_pieces.append(piece)

    metrics = pd.concat(met_pieces)
    print(metrics.columns)
    metrics = metrics.reset_index()
    metrics.rename(columns={'index': 'cell_id'}, inplace=True)

    metrics = metrics.loc[:, ~metrics.columns.duplicated()]

    # add columns for sample and library id
    metrics['library_id'] = metrics['cell_id'].apply(lambda x: x.split('-')[1])
    metrics['sample_id'] = metrics['cell_id'].apply(lambda x: x.split('-')[0])

    # Remap experimental conditions and filter
    conditions = {
        'A': 'G1',
        'A-BSA': 'G1',
        'A-NCC': 'G1',
        'B': 'S',
        'B-NCC': 'S',
        'C': 'G2',
        'C-NCC': 'G2',
        'G1': 'G1',
        'G2': 'G2',
        'S': 'S',
        'D': 'D',
    }

    conditions = pd.Series(conditions)
    conditions.index.name = 'experimental_condition'
    conditions.name = 'cell_cycle_state'
    conditions = conditions.reset_index()
    metrics = metrics.drop_duplicates()
    print(conditions.head())
    print(metrics.columns)
    print(metrics.head())
    metrics = metrics.merge(conditions, how='left')
    metrics = metrics.drop_duplicates()

    # merge the hmmcopy reads and annotation data
    cn = cn.merge(metrics, on='cell_id')
    cn = cn.loc[:, ~cn.columns.duplicated()]
    cn.drop_duplicates(inplace=True)

    # filter based on number of hmmcopy mapped reads
    cn = cn.query('total_mapped_reads_hmmcopy > 500000')

    # filter based on position (blacklisted loci have gc<0)
    cn = cn.query('gc > 0')

    # remove Y chromosome since we're dealing with female cell lines
    cn = cn.query('chr != "Y"')

    # filter out control cells that don't contain one of the sample_ids in their cell_id
    sample_ids = list(argv.samples)
    sample_ids.append(argv.dataset)
    # catch edge cases where cell_ids don't exactly match the sample or dataset id
    if argv.dataset in ['SA906a', 'SA906b']:
        sample_ids.append('SA906')
    elif argv.dataset == 'SA1292':
        sample_ids.append('AT135')
    elif 'SA609' in argv.dataset:
        sample_ids.append('SA609')
    elif 'SA1035' in argv.dataset:
        sample_ids.append('SA1035')
    elif 'SA535' in argv.dataset:
        sample_ids.append('SA535')
    print('sample_ids in list', sample_ids)
    print('sample_ids in cn', cn['sample_id'].unique())

    # remove cells that don't contain one of the sample_ids in their cell_id
    cn = cn[cn['cell_id'].str.contains('|'.join(sample_ids))]

    # filter out cells based on quality, experimental condition and contamination    
    cn = cn[
        (cn['cell_call'].isin(['C1', 'C2']))
    ]
    for rm_cond in ['gDNA', 'GM', 'NCC', 'NTC']:
        mask = ~cn['experimental_condition'].str.contains(rm_cond)
        cn = cn[mask]

    # convert is_contaminated to boolean before filtering
    cn['is_contaminated'] = cn['is_contaminated'].astype('bool')
    cn = cn[~cn['is_contaminated']]

    print("number of cells before filtering mouse bins", len(cn['cell_id'].unique()))
    cn = cn[cn['fastqscreen_mm10'] < cn['fastqscreen_grch37']]
    print("number of cells after filtering mouse bins", len(cn['cell_id'].unique()))

    # drop cell cycle state column because the mappings are inaccurate.. we only care about filtering on experimental condition
    cn.drop(columns=['cell_cycle_state'], inplace=True)

    # create new column for reads per million (normalize each cell's total read count to be the same)
    cn = compute_reads_per_million(cn, reads_col='reads', rpm_col='rpm')

    # save the output
    cn.to_csv(argv.output, sep='\t', index=False)
