from argparse import ArgumentParser
import numpy as np
import pandas as pd


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_input', help='input long-form copy number dataframe')
    p.add_argument('time_legend', help='CSV file containing timepoints for each library')
    p.add_argument('cn_out', help='output tsv that is same as cn_input with timepoint info added')

    return p.parse_args()


def main():
    argv = get_args()
    cn = pd.read_csv(argv.cn_input, sep='\t')
    time_df = pd.read_csv(argv.time_legend)

    # get rid of useless columns
    time_df = time_df[['library_id', 'label', 'datasetname', 'timepoint', 'sample_id']]

    # make sure chromosome column is set to the appropriate dtype
    cn['chr'] = cn['chr'].astype(str)
    cn['chr'] = cn['chr'].astype('category')

    # create library_id columns using cell_ids in the CN dataframes
    if 'library_id' not in cn.columns:
        cn['library_id'] = cn['cell_id'].apply(lambda x: x.split('-')[1])

    # merge time info with cn now using library_id as common column
    cn_out = pd.merge(cn, time_df)

    cn_out.to_csv(argv.cn_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
