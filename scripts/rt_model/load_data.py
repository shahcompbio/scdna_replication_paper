import pandas as pd
import numpy as np
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-t', '--table', type=str, help='table of input paths and metadata for each sample')
    p.add_argument('-s', '--samples', type=str, nargs='+', help='list of samples to include')
    p.add_argument('-cp', '--count_paths', type=str, nargs='+', help='cell cycle clone counts for across all samples (split across multiple files)')
    p.add_argument('-sr', '--sample_rt', type=str, help='matrix of RT profiles for each sample')
    p.add_argument('-cr', '--clone_rt', type=str, help='matrix of RT samples for each clone')
    p.add_argument('-sc', '--sample_cn', type=str, help='matrix of copy number profiles for each sample')
    p.add_argument('-cc', '--clone_cn', type=str, help='matrix of copy number profiles for each clone')
    p.add_argument('-f', '--features', type=str, help='clone feature metadata (e.g. signature, ploidy, cell type)')

    return p.parse_args()


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


def convert_to_one_hot(df, col):
    """ Given a dataframe and a column name, convert the column to one-hot encoding and return the new dataframe and the names of the new columns. """
    onehot = pd.get_dummies(df[col])
    onehot.columns = [col + '_' + t for t in onehot.columns]
    return pd.concat([df, onehot], axis=1).drop(col, axis=1), onehot.columns


def main():
    argv = get_args()

    # load the number of cells in each clone
    counts = load_counts(argv)

    # set a threshold for the number of S-phase cells in each clone
    num_cell_threshold = 10

    # load the table of input paths and metadata for each sample
    table = pd.read_csv(argv.table, sep='\t')

    # initialize empty dataframes for the sample and clone RT profiles
    sample_df = pd.DataFrame()
    clone_df = pd.DataFrame()
    sample_cn = pd.DataFrame()
    clone_cn = pd.DataFrame()

    # loop through each sample and load the RT and CN profiles
    for d in argv.samples:
        print('dataset:', d)
        temp_rt_path = table.query("dataset=='{}'".format(d))['rt_path'].values[0]
        temp_cn_path = table.query("dataset=='{}'".format(d))['cn_path'].values[0]
        if temp_rt_path.endswith('tsv'):
            temp_df = pd.read_csv(temp_rt_path, sep='\t')
            temp_cn = pd.read_csv(temp_cn_path, sep='\t')
        elif temp_rt_path.endswith('csv.gz'):
            temp_df = pd.read_csv(temp_rt_path)
            temp_cn = pd.read_csv(temp_cn_path)
        # drop all columns that contain '_hours' in the name
        temp_df = temp_df.loc[:, ~temp_df.columns.str.contains('_hours')]
        # drop all columns that end with '_U' or '_T' as we don't care about stratifying by treatment status
        # this is necessary for fitness samples as they have treated and untreated RT pseudobulk profiles
        temp_df = temp_df.loc[:, ~temp_df.columns.str.endswith('_U')]
        temp_df = temp_df.loc[:, ~temp_df.columns.str.endswith('_T')]
        # add 'end' as a column
        temp_df['end'] = temp_df['start'] + 500000 - 1
        temp_cn['end'] = temp_cn['start'] + 500000 - 1
        # set 'chr' and 'start' as index
        temp_df = temp_df.set_index(['chr', 'start', 'end'])
        # add the dataset name as the prefix to all column names
        temp_df = temp_df.add_prefix(d + '_')

        # subset to just the columns that contain 'clone' in the name
        temp_clone_df = temp_df.loc[:, temp_df.columns.str.contains('clone')]
        # subset to just the columns that don't contain 'clone' in the name
        temp_sample_df = temp_df.loc[:, ~temp_df.columns.str.contains('clone')]

        temp_cn = temp_cn.set_index(['chr', 'start', 'end'])
        temp_clone_cn = temp_cn.loc[:, temp_cn.columns.str.contains('clone')]
        temp_clone_cn = temp_clone_cn.add_prefix(d + '_')
        temp_sample_cn = temp_cn.loc[:, temp_cn.columns.str.contains('sample')]

        # loop through each clone and check to see if they have enough S-phase cells
        good_clones = []
        good_clone_cn_cols = []
        for col in temp_clone_df.columns:
            print('col:', col)
            c = col.split('_')[2].replace('clone', '')
            print('clone:', c)
            clone_cn_col = '{d}_clone_{c}'.format(d=d, c=c)
            try:
                n = int(counts.query("dataset=='{}'".format(d)).query("clone_id=='{}'".format(c))['num_cells_s'].values[0])
            except:
                n = 0  # if a clone is missing from 'counts', it means there are no S-phase cells
            if n > num_cell_threshold:
                good_clones.append(col)
                good_clone_cn_cols.append(clone_cn_col)
        # subset the clone df to just the good clones
        temp_clone_df = temp_clone_df[good_clones]
        temp_clone_cn = temp_clone_cn[good_clone_cn_cols]

        # reset the index so that chr, start are columns again
        temp_sample_df = temp_sample_df.reset_index()
        temp_clone_df = temp_clone_df.reset_index()
        temp_sample_cn = temp_sample_cn.reset_index()
        temp_clone_cn = temp_clone_cn.reset_index()

        # if df is empty, set it to temp_df, otherwise merge temp_df with df
        # start with the sample df
        if sample_df.empty:
            sample_df = temp_sample_df
        else:
            sample_df = sample_df.merge(temp_sample_df)

        # do the same for the clone df
        if clone_df.empty:
            clone_df = temp_clone_df
        else:
            clone_df = clone_df.merge(temp_clone_df)
        
        # do the same for the sample cn
        if sample_cn.empty:
            sample_cn = temp_sample_cn
        else:
            sample_cn = sample_cn.merge(temp_sample_cn)
        
        # do the same for the clone cn
        if clone_cn.empty:
            clone_cn = temp_clone_cn
        else:
            clone_cn = clone_cn.merge(temp_clone_cn)
    
    # rename all columns with 'dataset' prefix to 'sample' prefix in sample_cn
    sample_cn = sample_cn.rename(columns={c: c.replace('dataset', 'sample') for c in sample_cn.columns})

    # add an empty column in counts for 'ploidy'
    counts['ploidy'] = np.nan
    # iterate through each row of counts
    for i, row in counts.iterrows():
        c, d = row['clone_id'], row['dataset']
        col = '{d}_clone_{c}'.format(d=d, c=c)
        if col in clone_cn.columns:
            # get the ploidy value from the clone_cn df
            ploidy = np.median(clone_cn['{d}_clone_{c}'.format(d=d, c=c)].values)
            # set the ploidy value in counts
            counts.loc[i, 'ploidy'] = ploidy

    # merge counts and sample mapping
    features = pd.merge(counts, table.drop(columns=['rt_path', 'cn_path']), on='dataset')

    # convert NaN signatures into string 'N/A' so they don't get filtered out
    features['signature'] = features['signature'].fillna('N/A')

    # drop rows with ploidy==NaN in features but leave 
    print(features.shape)
    features = features.dropna()
    print(features.shape)

    # set chr, start, end as indices for sample_cn, sample_df, clone_cn, clone_df
    clone_cn.set_index(['chr', 'start', 'end'], inplace=True)
    clone_df.set_index(['chr', 'start', 'end'], inplace=True)
    sample_cn.set_index(['chr', 'start', 'end'], inplace=True)
    sample_df.set_index(['chr', 'start', 'end'], inplace=True)

    # remove the genomic loci that are found in one clone dataframe but not the other
    clone_cn_cols = clone_cn.columns
    clone_df_cols = clone_df.columns
    temp_merged_clone = pd.merge(clone_cn.reset_index(), clone_df.reset_index())
    temp_merged_clone.set_index(['chr', 'start', 'end'], inplace=True)
    clone_cn = temp_merged_clone[clone_cn_cols]
    clone_df = temp_merged_clone[clone_df_cols]
    print(clone_cn.shape, clone_df.shape)

    # create a list of categorical columns to convert to one-hot encoding
    categorical_columns = ['type', 'signature']
    # create a list of all feature columns
    # initialize with the non-catagorical features (ploidy, etc.)
    feature_cols = ['ploidy']
    clone_features = features.copy()
    # convert each column to one-hot encoding
    for col in categorical_columns:
        clone_features, new_cols = convert_to_one_hot(clone_features, col)
        feature_cols.extend(new_cols)
    
    print('feature_cols:', feature_cols)

    # remove signature_Line from feature_cols and clone_features since it is just a dummy variable for hTERT samples with no signature
    feature_cols.remove('signature_N/A')
    clone_features = clone_features.drop(['signature_N/A'], axis=1)

    # print the shape of files that I'll be saving
    print('clone_df:', clone_df.shape)
    print('sample_df:', sample_df.shape)
    print('clone_cn:', clone_cn.shape)
    print('sample_cn:', sample_cn.shape)
    print('clone_features:', clone_features.shape)

    # save ouput files
    clone_df.to_csv(argv.clone_rt)
    sample_df.to_csv(argv.sample_rt)
    clone_cn.to_csv(argv.clone_cn)
    sample_cn.to_csv(argv.sample_cn)
    clone_features.to_csv(argv.features, index=False)


if __name__=='__main__':
    main()