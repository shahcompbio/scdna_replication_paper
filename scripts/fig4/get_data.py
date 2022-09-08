import logging
import pandas as pd
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('--cn_data_urls', type=str, nargs='+', help='list of reads data from different libraries')
    p.add_argument('--metrics_data_urls', type=str, nargs='+', help='list of metrics data from different libraries')
    p.add_argument('--align_metrics_data_urls', type=str, nargs='+', help='list of align metrics data from different libraries')
    p.add_argument('--curated_labels', nargs='?', default=None, help='file containing curated cell cycle state labels (optional)')
    p.add_argument('--cn_out', type=str, help='output for all cn data')
    p.add_argument('--metrics_out', type=str, help='output for all metrics data')
    p.add_argument('--align_metrics_out', type=str, help='output for all align metrics data')

    return p.parse_args()


def get_data(argv, sas='',):
    cn_data = []
    for cn_data_url in argv.cn_data_urls:
        url = cn_data_url + sas
        logging.info(f'reading from {url}')
        cn_data.append(pd.read_csv(url, compression='gzip'))
    cn_data = pd.concat(cn_data, sort=True, ignore_index=True)

    metrics_data = []
    for metrics_data_url in argv.metrics_data_urls:
        url = metrics_data_url + sas
        logging.info(f'reading from {url}')
        metrics_data.append(pd.read_csv(url, compression='gzip'))
    metrics_data = pd.concat(metrics_data, sort=True, ignore_index=True)

    align_metrics_data = []
    for align_metrics_data_url in argv.align_metrics_data_urls:
        url = align_metrics_data_url + sas
        logging.info(f'reading from {url}')
        align_metrics_data.append(pd.read_csv(url, compression='gzip'))
    align_metrics_data = pd.concat(align_metrics_data, sort=True, ignore_index=True)

    for data in (cn_data, metrics_data, align_metrics_data):
        data['sample_id'] = [a.split('-')[0] for a in data['cell_id']]
        data['library_id'] = [a.split('-')[1] for a in data['cell_id']]

    # Fix total_mapped_reads_hmmcopy column
    fix_read_count = metrics_data['total_mapped_reads_hmmcopy'].isnull()
    metrics_data.loc[fix_read_count, 'total_mapped_reads_hmmcopy'] = (
        metrics_data.loc[fix_read_count, 'total_mapped_reads'])

    metrics_data = metrics_data.query('total_mapped_reads_hmmcopy > 500000')

    logging.info('library sizes:\n{}'.format(metrics_data.groupby('library_id').size()))

    cn_data = cn_data.merge(metrics_data[['cell_id']].drop_duplicates())
    cn_data = cn_data.merge(metrics_data[['cell_id', 'experimental_condition', 'log_likelihood', 'MBRSM_dispersion',
                                        'MBRSI_dispersion_non_integerness', 'mean_state_mads', 'mad_hmmcopy']])

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

    metrics_data = metrics_data.merge(conditions)
    cn_data = cn_data.merge(conditions)

    metrics_data = metrics_data.query('cell_cycle_state != "D"')
    cn_data = cn_data.query('cell_cycle_state != "D"')

    # have hand-curated labels replace the flow-determined cell cycle state
    if argv.curated_labels is not None:
        curated_df = pd.read_csv(argv.curated_labels)

        # merge new_cell_cycle_state into metrics data
        metrics_data = pd.merge(metrics_data, curated_df, how='left', on=['cell_id', 'cell_cycle_state'])

        # copy flow-determined state for cells without new_cell_cycle_state
        for key, row in metrics_data.iterrows():
            if pd.isnull(row['new_cell_cycle_state']):
                metrics_data.loc[key, 'new_cell_cycle_state'] = metrics_data.loc[key, 'cell_cycle_state']
                metrics_data.loc[key, 'checked'] = False
                metrics_data.loc[key, 'corrected'] = False

        # merge assignments into cn_data
        full_curated_df = metrics_data[['cell_id', 'new_cell_cycle_state', 'checked', 'corrected']]
        cn_data = pd.merge(cn_data, full_curated_df, how='left', on='cell_id')

        # assign new_cell_cycle_state to cell_cycle_state
        cn_data.cell_cycle_state = cn_data.new_cell_cycle_state
        metrics_data.cell_cycle_state = metrics_data.new_cell_cycle_state
        cn_data.drop(columns=['new_cell_cycle_state'], inplace=True)
        metrics_data.drop(columns=['new_cell_cycle_state'], inplace=True)

    return cn_data, metrics_data, align_metrics_data


if __name__ == '__main__':
    argv = get_args()
    cn_data, metrics_data, align_metrics_data = get_data(argv)

    cn_data.to_csv(argv.cn_out, sep='\t', index=False)
    metrics_data.to_csv(argv.metrics_out, sep='\t', index=False)
    align_metrics_data.to_csv(argv.align_metrics_out, sep='\t', index=False)
