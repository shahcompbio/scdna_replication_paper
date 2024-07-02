from argparse import ArgumentParser
import pandas as pd


def get_args():
    p = ArgumentParser()

    p.add_argument('--input', type=str, help='csv file of all the signatures paths in isabl')
    p.add_argument('--dataset', type=str, help='name of this dataset -- used for filtering out control cells')
    p.add_argument('--output', type=str, help='output tsv file that combines all cells from same dataset into one file')

    return p.parse_args()


def compute_reads_per_million(cn, reads_col='reads', rpm_col='rpm'):
    for cell_id, cell_cn in cn.groupby('cell_id'):
        x = cell_cn[reads_col].values
        temp_rpm = x * 1e6 / sum(x)
        cn.loc[cell_cn.index, rpm_col] = temp_rpm
    return cn


def dataset_cn_files(isabl_table, dataset):
    files = isabl_table.loc[
        isabl_table['isabl_patient_id']==dataset].loc[
        isabl_table['result_type']=='reads']['result_filepath'].values
    return list(files)


def dataset_metric_files(isabl_table, dataset):
    files = isabl_table.loc[
        isabl_table['isabl_patient_id']==dataset].loc[
        isabl_table['result_type']=='metrics']['result_filepath'].values
    return list(files)


def dataset_sample_ids(isabl_table, dataset):
    sample_ids = isabl_table.loc[isabl_table['isabl_patient_id']==dataset]['isabl_sample_id'].unique()
    return list(sample_ids)


if __name__ == '__main__':
    argv = get_args()

    isabl_table = pd.read_csv(argv.input)

    hmm_files = dataset_cn_files(isabl_table, argv.dataset)
    met_files = dataset_metric_files(isabl_table, argv.dataset)
    sample_ids = dataset_sample_ids(isabl_table, argv.dataset)

    cn_pieces = []

    for f in hmm_files:
        piece = pd.read_csv(f)
        piece = piece[['chr', 'start', 'end', 'cell_id', 'reads', 'gc', 'map', 'copy', 'state']]
        cn_pieces.append(piece)

    cn = pd.concat(cn_pieces, ignore_index=True)

    met_pieces = []
    for f, isabl_sample_id in zip(met_files, sample_ids):
        piece = pd.read_csv(f)
        piece = piece[[
            'total_mapped_reads_hmmcopy', 'is_s_phase', 'is_s_phase_prob', 'multiplier',
            'quality', 'is_contaminated', 'fastqscreen_human', 'fastqscreen_mouse', 'coverage_depth',
            'breakpoints', 'sample_id', 'library_id', 'cell_id'
        ]]
        piece['isabl_sample_id'] = isabl_sample_id
        met_pieces.append(piece)

    metrics = pd.concat(met_pieces, ignore_index=True)
    print(metrics.columns)

    # filter out duplicate columns
    metrics = metrics.loc[:, ~metrics.columns.duplicated()]

    # list of control cell sample_id values that need to be filtered out
    control_sample_ids = [
        'DLPNegative',
        'DLPGm'
    ]
    # filter out control cells
    metrics = metrics[~metrics['sample_id'].isin(control_sample_ids)]

    # merge metrics and cn dataframes
    cn = cn.merge(metrics, on='cell_id')
    cn = cn.loc[:, ~cn.columns.duplicated()]
    cn.drop_duplicates(inplace=True)

    # add the isabl_sample_id as the prefix to the cell_id
    cn['cell_id'] = cn['isabl_sample_id'] + '-' + cn['cell_id']

    # add a column for the isabl_patient_id
    cn['isabl_patient_id'] = argv.dataset

    # filter based on number of hmmcopy mapped reads
    cn = cn.query('total_mapped_reads_hmmcopy > 500000')

    # filter based on position (blacklisted loci have gc<0)
    cn = cn.query('gc > 0')

    # remove Y chromosome since we're dealing with female patients
    cn = cn.query('chr != "Y"')

    # convert is_contaminated to boolean and filter out contaminated cells
    cn['is_contaminated'] = cn['is_contaminated'].astype('bool')
    cn = cn[~cn['is_contaminated']]

    # filter out cells with more mouse than human reads
    print("shape before filtering mouse bins", cn.shape)
    cn = cn[cn['fastqscreen_mouse'] < cn['fastqscreen_human']]
    print("shape after filtering mouse bins", cn.shape)

    # create new column for reads per million (normalize each cell's total read count to be the same)
    cn = compute_reads_per_million(cn, reads_col='reads', rpm_col='rpm')

    cn.to_csv(argv.output, index=False)
