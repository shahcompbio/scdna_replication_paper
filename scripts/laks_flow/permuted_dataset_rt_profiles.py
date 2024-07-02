from argparse import ArgumentParser
import pandas as pd


def get_args():
    p = ArgumentParser()

    p.add_argument('-rr', '--rt_ref', type=str, help='rt pseudobulks when no labels are permuted')
    p.add_argument('-rp', '--rt_perm', type=str, nargs='+', help='rt pseudobulks for all the permuted datasets')
    p.add_argument('-d', '--datasets', type=str, nargs='+', help='dataset ids')
    p.add_argument('-rt', '--rt_table', help='table containing all the rt pseudobulk profiles for all permuted datasets')

    return p.parse_args()


def load_data(argv, i, d):
    rt = pd.read_csv(argv.rt_perm[i], sep='\t')
    rt['dataset'] = d
    return rt


def main():
    argv = get_args()

    # load rt bulks from dataset with no permuted labels (use as a reference)
    ref_rt = pd.read_csv(argv.rt_ref, sep='\t')

    # build table that matches config file for each permuted dataset
    legend_df = pd.DataFrame({
        'dataset': argv.datasets,
    })

    # load in all the permuted dataset bulk rt profiles
    all_rts = []
    i = 0
    for dataset, row in legend_df.groupby('dataset'):
        temp_rt = load_data(argv, i, dataset)
        all_rts.append(temp_rt)
        i += 1
    all_rts = pd.concat(all_rts, ignore_index=True)

    # merge each permuted dataset's relevant RT columns into rt_wide
    rt_wide = ref_rt.copy()
    bulk_rt_cols = ['rt_merged_T47D', 'rt_merged_GM18507']

    for dataset, temp_rt in all_rts.groupby('dataset'):
        # rename columns in temp_rt to reflect the current dataset
        t_col = 'T47D_{}'.format(dataset)
        gm_col = 'GM18507_{}'.format(dataset)
        bulk_rt_cols.extend([t_col, gm_col])
        temp_rt = temp_rt.rename(columns={
            'rt_T47D': t_col,
            'rt_GM18507': gm_col,
        })
        
        # merge the columns using loci
        temp_rt = temp_rt[['chr', 'start', t_col, gm_col]]
        rt_wide = pd.merge(rt_wide, temp_rt)
    
    # save the table used to compute the correlations
    rt_wide.to_csv(argv.rt_table, sep='\t', index=False)


if __name__ == '__main__':
    main()
