import pandas as pd
import numpy as np
from scdna_replication_tools.compute_pseudobulk_rt_profiles import compute_pseudobulk_rt_profiles
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-t', '--times', type=str, help='path to csv file with time points of all the fitness libraries (including SA039, SA906a, SA906b)')
    p.add_argument('-s', '--scRT', type=str, nargs='+', help='full length PERT dataframes from SA039, SA906a, SA906b')
    p.add_argument('-p', '--profiles', type=str, help='subsampled pseudobulk RT profiles for each sample, clone, time triplet')
    p.add_argument('-f', '--features', type=str, help='features that correspond to the pseudobulk RT profiles')

    return p.parse_args()


def subsampled_rt_profiles(scrt_df, cell_counts, num_cells=25, num_samples=100, cell_thresh=50):
    ''' 
    Given a dataframe of PERT results and a set of [datasetname, clone_id, timepoint] triplets,
    subsample num_cells from each triplet and compute the mean RT profile for each triplet. 
    '''
    # subset cell counts to just the triplets with sufficeint cells
    cell_counts = cell_counts[cell_counts['num_s_cells'] >= cell_thresh]

    # loop through all the triplets
    rt_profiles = pd.DataFrame()
    features_df = []
    for i, row in cell_counts.iterrows():
        d = row['datasetname']
        c = row['clone_id']
        t = row['timepoint']
        # subset the scRT dataframe to just the cells from this triplet
        df = scrt_df[(scrt_df['datasetname'] == d) & (scrt_df['clone_id'] == c) & (scrt_df['timepoint'] == t)]
        temp_cell_ids = df['cell_id'].unique()
        for j in range(num_samples):
            np.random.seed(j)
            # randomly sample num_cells from the cells in this triplet
            cell_ids = np.random.choice(temp_cell_ids, size=num_cells, replace=False)
            # subset the scRT dataframe to just the sampled cells
            temp_df = df[df['cell_id'].isin(cell_ids)]
            # compute the pseudobulk RT profile for the sampled cells
            temp_rt = compute_pseudobulk_rt_profiles(temp_df, 'model_rep_state', clone_col=None, output_col='pseudobulk')
            # subset to just chr, start, and pseudobulk_model_rep_state
            temp_rt = temp_rt[['chr', 'start', 'pseudobulk_model_rep_state']]
            # rename the pseudobulk_model_rep_state column to include the datasetname, clone_id, timepoint and the sample number
            temp_rt.columns = ['chr', 'start', '{}_{}_{}_{}'.format(d, c, t, j)]
            # merge the pseudobulk RT profile for this sample with the others
            if rt_profiles.empty:
                rt_profiles = temp_rt
            else:
                rt_profiles = pd.merge(rt_profiles, temp_rt, on=['chr', 'start'])
            # store the features for this sample
            temp_features = pd.DataFrame({
                'datasetname': [d], 'clone_id': [c], 'timepoint': [t], 'sample': [j], 'subsampled_cells': [num_cells], 'total_cells': [row['num_s_cells']],
            })
            features_df.append(temp_features)
    
    # merge the features for all the samples
    features_df = pd.concat(features_df, ignore_index=True)

    return rt_profiles, features_df


def main():
    argv = get_args()

    # load data that maps library id to timepoints
    time_df = pd.read_csv(argv.times)
    datasetnames = ['SA039U', 'SA906a', 'SA906b']
    time_df = time_df[time_df['datasetname'].isin(datasetnames)]
    # rename datasetname with 'SA039U' to 'SA039'
    time_df.loc[time_df['datasetname'] == 'SA039U', 'datasetname'] = 'SA039'
    # subset to the necessary columns
    time_df = time_df[['datasetname', 'timepoint', 'library_id']]

    # load the scRT results for each dataset
    scrt_df = []
    for path in argv.scRT:
        d = path.split('/')[-2]
        df = pd.read_csv(path, sep='\t')
        print('done loading dataset', d)
        df['datasetname'] = d
        scrt_df.append(df)
    scrt_df = pd.concat(scrt_df, ignore_index=True)

    # merge the timepoints with the scRT results using the library_id
    scrt_df = pd.merge(scrt_df, time_df, on=['datasetname', 'library_id'])

    # create a new dataset+clone+timepoint column to use for grouping
    scrt_df['dct'] = scrt_df['datasetname'] + '_' + scrt_df['clone_id'].astype(str) + '_' + scrt_df['timepoint'].astype(str)

    # compute the number of S-phase cells at each timepoint & time
    cell_counts = scrt_df[['datasetname', 'clone_id', 'timepoint', 'dct', 'cell_id']].drop_duplicates()[['datasetname', 'clone_id', 'timepoint', 'dct']].value_counts().to_frame().reset_index().rename(columns={0: 'num_s_cells'})

    # compute subsampled pseudobulk RT profiles for each dataset+clone+timepoint group
    rt_profiles, features_df = subsampled_rt_profiles(scrt_df, cell_counts)

    # save the rt_profiles and features_df
    rt_profiles.to_csv(argv.profiles, index=False)
    features_df.to_csv(argv.features, index=False)


if __name__ == '__main__':
    main()
