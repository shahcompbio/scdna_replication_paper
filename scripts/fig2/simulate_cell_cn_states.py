import numpy as np
import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-gr', '--ref_data', help='gc and rt data at different bin sizes')
    p.add_argument('-gm', '--gc_profile_500kb', help='gc and map data at 500kb without blacklisted regions')
    p.add_argument('-nS', '--num_cells_S', type=int, help='number of S-phase cells')
    p.add_argument('-nG', '--num_cells_G', type=int, help='number of G1/2-phase cells')
    p.add_argument('-bs', '--bin_size', type=int, help='bin size to use when simulating data (must be in ref_data)')
    p.add_argument('-c', '--clones', type=str, nargs='+', help='list of unique clone ids')
    p.add_argument('-cp', '--clone_probs', type=float, nargs='+', help='probability of observing each clone')
    p.add_argument('-s', '--states', type=int, nargs='+', help='list of unique copy number states')
    p.add_argument('-sp', '--state_probs', type=float, nargs='+', help='probability of each copy number state occurring')
    p.add_argument('-cna', '--cell_CNA_prob', type=float, help='probability of each chromosome having a cell-specific CNA')
    p.add_argument('-rt', '--rt_col', type=str, help='column in ref_data containing appropriate replication timing values')
    p.add_argument('-so', '--s_out', help='simulated S-phase cells')
    p.add_argument('-go', '--g_out', help='simulated G1/2-phase cells')

    return p.parse_args()


def filter_small_bin_df(large_bin_df, small_bin_df):
    """ Remove any small bins that don't lie inside a large bin. """
    filtered_small_df = []
    for i, row in large_bin_df.iterrows():
        chrom = row.chr
        start = row.start - 1
        end = row.end
        temp = small_bin_df.loc[(small_bin_df['chr']==chrom) & (small_bin_df['start']>=start) & (small_bin_df['end']<=end)]
        filtered_small_df.append(temp)
    filtered_small_df = pd.concat(filtered_small_df, ignore_index=True)
    return filtered_small_df


def simulate_clone_profiles(ref_data, clones, states, state_probs, clone_cn_column='clone_cn_state'):
    clone_cn = []
    for clone_id in clones:
        clone_df = ref_data.copy()
        clone_df['clone_id'] = clone_id
        clone_df.reset_index(inplace=True, drop=True)
        for chrom_id, chunk in clone_df.groupby('chr'):
            cn = np.random.choice(states, p=state_probs)
            clone_df.loc[chunk.index, clone_cn_column] = cn
        clone_cn.append(clone_df)
    clone_cn = pd.concat(clone_cn, ignore_index=True)
    return clone_cn


def simulate_cell_profiles(clone_cn, num_cells, clones, clone_probs, states, state_probs, cell_CNA_prob, rt_col='mcf7rt',
                           cell_id_prefix='cell_S', clone_cn_column='clone_cn_state', cell_cn_column='true_G1_state'):
    cell_cn = []
    for i in range(num_cells):
        # draw the clone this cell is from
        clone_id = np.random.choice(clones, p=clone_probs)
        temp_clone_df = clone_cn.loc[clone_cn['clone_id']==clone_id]
        # create a cell-specific copy of the clone profile
        temp_cell_df = temp_clone_df.copy()
        temp_cell_df['cell_id'] = '{}_{}'.format(cell_id_prefix, i)
        temp_cell_df[cell_cn_column] = temp_cell_df[clone_cn_column]
        temp_cell_df.reset_index(inplace=True, drop=True)
        # loop through each chromosome and add some cell-specific CNAs
        for chrom_id, chunk in temp_cell_df.groupby('chr'):
            change_cn = np.random.choice([0, 1], p=[1-cell_CNA_prob, cell_CNA_prob])
            current_cn = chunk[cell_cn_column].values[0]
            # find a new CN value if true
            if change_cn == 1:
                keep_going = True
                while keep_going:
                    cn = np.random.choice(states, p=state_probs)
                    if cn != current_cn:
                        keep_going = False
                temp_cell_df.loc[chunk.index, cell_cn_column] = cn

        # add cell to list
        cell_cn.append(temp_cell_df)
    cell_cn = pd.concat(cell_cn, ignore_index=True)

    return cell_cn


def main():
    argv = get_args()

    ref_data = pd.read_csv(argv.ref_data, dtype={'Chromosome': str})
    ref_data = ref_data.dropna().reset_index(drop=True)
    ref_data.rename(columns={'Chromosome': 'chr', 'Start': 'start', 'End': 'end'}, inplace=True)
    ref_data['chr'] = ref_data['chr'].astype('category')

    ref_data = ref_data.loc[ref_data['bin_size']==int(argv.bin_size)]
    ref_data.reset_index(inplace=True, drop=True)
    print('ref_data.head()\n', ref_data.head())

    gc_profile_500kb = pd.read_csv(argv.gc_profile_500kb)
    ref_data = filter_small_bin_df(gc_profile_500kb, ref_data)
    print('ref_data.head()\n', ref_data.head())

    assert len(argv.clones) == len(argv.clone_probs)
    assert len(argv.states) == len(argv.state_probs)

    # simulate consensus clone profiles
    clone_cn = simulate_clone_profiles(ref_data, argv.clones, argv.states, argv.state_probs)

    print('simulating S-phase cells...')
    s_cells = simulate_cell_profiles(
        clone_cn, argv.num_cells_S, argv.clones, argv.clone_probs, argv.states, argv.state_probs, argv.cell_CNA_prob,
        cell_id_prefix='cell_S', rt_col=argv.rt_col
    )
    print('s_cells.head()\n', s_cells.head())

    print('simulating G1-phase cells...')
    g_cells = simulate_cell_profiles(
        clone_cn, argv.num_cells_G, argv.clones, argv.clone_probs, argv.states, argv.state_probs, argv.cell_CNA_prob,
        cell_id_prefix='cell_G'
    )
    print('g_cells.head()\n', g_cells.head())

    s_cells.to_csv(argv.s_out, sep='\t', index=False)
    g_cells.to_csv(argv.g_out, sep='\t', index=False)



if __name__ == '__main__':
    main()