from argparse import ArgumentParser
import pandas as pd
from scdna_replication_tools.compute_ccc_features import compute_ccc_features


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', help='input long-form copy number dataframe for cells in all cell cycle phases')
    p.add_argument('cn_out', help='same as cn_input with classifier features added')


    return p.parse_args()


def main():
    argv = get_args()

    # load in data
    cn = pd.read_csv(argv.cn_input)

    # treat each library as a unique clone when computing the cell cycle classifier features
    cn_temp, cell_features = compute_ccc_features(
        cn, cell_col='cell_id', rpm_col='rpm', clone_col='clone_id', madn_col='madn',
        lrs_col='lrs', num_reads_col='total_mapped_reads_hmmcopy', bk_col='breakpoints'
    )

    # merge the cell level features with the cn input
    # we don't want to use cn_temp as the output since its
    # shape might be different from the cn_input which has already been filtered
    cn_out = pd.merge(cn, cell_features)

    # save output files
    cn_out.to_csv(argv.cn_out, index=False)


if __name__ == '__main__':
    main()
