from argparse import ArgumentParser
import numpy as np
import pandas as pd


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s_input', help='input long-form copy number dataframe for s-phase cells')
    p.add_argument('cn_g_input', help='input long-form copy number dataframe for g1/2-phase cells')
    p.add_argument('time_legend', help='CSV file containing timepoints for each library')
    p.add_argument('cn_s_out', help='output tsv that is same as cn_s_input with timepoint info added')
    p.add_argument('cn_g_out', help='output tsv that is same as cn_g_input with timepoint info added')

    return p.parse_args()


def main():
    argv = get_args()
    cn_s = pd.read_csv(argv.cn_s_input, sep='\t')
    cn_g = pd.read_csv(argv.cn_g_input, sep='\t')
    time_df = pd.read_csv(argv.time_legend)

    # get rid of useless columns
    time_df = time_df[['library_id', 'label', 'timepoint']]

    # manually filter out low quality libraries
    bad_libraries = ['A96216A']
    time_df = time_df.loc[~time_df['library_id'].isin(bad_libraries)]

    # make sure chromosome column is set to the appropriate dtype
    cn_s['chr'] = cn_s['chr'].astype(str)
    cn_g['chr'] = cn_g['chr'].astype(str)

    # create library_id columns using cell_ids in the CN dataframes
    if 'library_id' not in cn_s.columns:
        cn_s['library_id'] = cn_s['cell_id'].apply(lambda x: x.split('-')[1])
    if 'library_id' not in cn_g.columns:
        cn_g['library_id'] = cn_g['cell_id'].apply(lambda x: x.split('-')[1])

    # merge time info with cn now using library_id as common column
    cn_s_out = pd.merge(cn_s, time_df)
    cn_g_out = pd.merge(cn_g, time_df)

    # filter out libraries that only have a few

    cn_s_out.to_csv(argv.cn_s_out, sep='\t', index=False)
    cn_g_out.to_csv(argv.cn_g_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
