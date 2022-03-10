import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('ref_data', help='gc and rt data at different bin sizes')
    p.add_argument('gc_profile_500kb', help='gc and map data at 500kb without blacklisted regions')
    p.add_argument('sigma1', type=float, help='noise of read depth profiles')
    p.add_argument('gc_slope', type=float, help='slope of linear GC bias')
    p.add_argument('gc_int', type=float, help='intercept of linear GC bias')
    p.add_argument('A', type=float, help='steepness of inflection point when drawing RT state')
    p.add_argument('B', type=float, help='offset for where inflection point occurs when drawing RT state')
    p.add_argument('num_reads', type=int, help='number of reads per cell')
    p.add_argument('num_cells_S', type=int, help='number of S-phase cells')
    p.add_argument('num_cells_G', type=int, help='number of G1/2-phase cells')
    p.add_argument('bin_size', type=int, help='bin size to use when simulating data (must be in ref_data)')
    p.add_argument('s_time_stdev', type=float, help='standard deviation for normal distribution of S-phase times/fractions (negative for uniform distribution)')
    p.add_argument('s_out', help='simulated S-phase cells')
    p.add_argument('g_out', help='simulated G1/2-phase cells')

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


def softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0)


def model(G1_state, gc, rt, s_time, sigma1, gc_slope, gc_int, num_reads, A=1, B=0, allele_specific=False):
    """Given the true CN state, GC, and RT values for a cell, come up with that observed read count profile."""
    # probablility of each allele being replicated follows logistic function of rt minus s_time
    # A controls the steepness; large A --> steeper inflection point --> deterministic replication
    # B determines where the inflection point occurs
    p_rep = 1 / (1 + np.exp(-A*(rt - s_time - B)))
    if allele_specific:
        # each allele can be replicated up to G1_state times
        replicated = np.random.binomial(G1_state, p_rep, rt.shape[0])
        true_CN = G1_state + replicated
    else:
        # assume that all alleles at each loci are either replicated or unreplicated
        replicated = np.random.binomial(1, p_rep, rt.shape[0])
        true_CN = G1_state * (1 + replicated)
    
    # add gc bias to the true CN
    # Is a simple linear model sufficient here?
    observed_CN = true_CN * ((gc * gc_slope) + gc_int)
    
    # add some random noise to the observed copy number
    reads_norm = softplus(np.random.normal(loc=observed_CN, scale=sigma1))
    
    # scale reads_norm and then draw true read count from multinomial distribution
    expected_reads_pval = reads_norm / sum(reads_norm)
    read_count = np.random.multinomial(num_reads, expected_reads_pval)
    
    return read_count, replicated, true_CN, observed_CN


def simulate_cell(temp_df, cell_id, s_time, num_reads, gc_slope, gc_int, sigma1, A, B):    
    read_count, replicated, true_CN, observed_CN = model(
        temp_df['true_G1_state'].values, temp_df['gc'].values, temp_df['mcf7rt'].values,
        s_time, sigma1, gc_slope, gc_int, num_reads, A=A, B=B
    )
    
    # store simulated cell values in temp_df
    temp_df['reads'] = read_count
    temp_df['true_rt_state'] = replicated
    temp_df['true_observed_CN'] = observed_CN
    temp_df['true_s_time'] = s_time
    temp_df['true_frac_rt'] = sum(replicated) / len(replicated)
    # exact number of reads will be a few off from input due to rounding errors --> use sum of output
    temp_df['total_mapped_reads_hmmcopy'] = sum(read_count)
    temp_df['gc_slope'] = gc_slope
    temp_df['gc_int'] = gc_int
    temp_df['sigma1'] = sigma1
    temp_df['A'] = A
    temp_df['B'] = B
    temp_df['cell_id'] = cell_id
    temp_df['rpm'] = (temp_df['reads'] / sum(temp_df['reads'].values)) * 1e6
    
    return temp_df


def simulate_diploid_s_cells(
    ref_data, num_cells, gc_slope, gc_int, sigma1, A, B, num_reads,
    rt_col='mcf7rt', cell_id_prefix='cell_S', s_time_stdev=10
):
    sim_s_df = []
    for i in range(num_cells):
        temp_df = ref_data.copy()
        cell_id = '{}_{}'.format(cell_id_prefix, i)

        min_rt_value = min(temp_df[rt_col].values)
        max_rt_value = max(temp_df[rt_col].values)
        
        if s_time_stdev <= 0:
            # draw from a normal distribution when stdev is not positive
            s_time = np.random.uniform(low=min_rt_value, high=max_rt_value)
        else:
            # draw from normal distribution and check too see if within range of min-max of rt_col
            repeat = True
            while repeat:
                s_time = np.random.normal(loc=np.mean(temp_df[rt_col].values), scale=s_time_stdev)
                if s_time < max_rt_value and s_time > min_rt_value:
                    repeat = False

        # assume a diploid cell
        temp_df['true_G1_state'] = 2
        # simulate gc bias and to get read counts
        temp_df = simulate_cell(temp_df, cell_id, s_time, num_reads, gc_slope, gc_int, sigma1, A, B)

        sim_s_df.append(temp_df)

    sim_s_df = pd.concat(sim_s_df, ignore_index=True)
    
    return sim_s_df


def simulate_diploid_g_cells(
    ref_data, num_cells, gc_slope, gc_int, sigma1, A, B, num_reads,
    rt_col='mcf7rt', cell_id_prefix='cell_S'
):
    sim_g_df = []
    for i in range(num_cells):
        temp_df = ref_data.copy()
        cell_id = '{}_{}'.format(cell_id_prefix, i)

        s_time = np.inf

        # assume a diploid cell
        temp_df['true_G1_state'] = 2
        # simulate gc bias and to get read counts
        temp_df = simulate_cell(temp_df, cell_id, s_time, num_reads, gc_slope, gc_int, sigma1, A, B)

        sim_g_df.append(temp_df)

    sim_g_df = pd.concat(sim_g_df, ignore_index=True)
    
    return sim_g_df


def aggregate_small_bin_df(large_bin_df, small_bin_df):
    """ Aggregate a small bin df into larger bins. """
    aggregated_df = []
    for i, row in large_bin_df.iterrows():
        chrom = row.chr
        start = row.start - 1
        end = row.end
        # extract small bins within this window
        temp = small_bin_df.loc[(small_bin_df['chr']==chrom) & (small_bin_df['start']>=start) & (small_bin_df['end']<=end)]
        
        # aggregate reads, copy, state, and replicated 
        reads = sum(temp['reads'])
        true_observed_CN = np.mean(temp['true_observed_CN'])
        state = np.median(temp['true_G1_state'])
        rt_state = np.median(temp['true_rt_state'])
        rt_norm = np.mean(temp['true_rt_state'])
        
        # per-cell values
        cell_id = temp.cell_id.values[0]
        true_s_time = temp.true_s_time.values[0]
        gc_slope = temp.gc_slope.values[0]
        gc_int = temp.gc_int.values[0]
        sigma1 = temp.sigma1.values[0]
        A = temp.A.values[0]
        B = temp.B.values[0]
        true_frac_rt = temp.true_frac_rt.values[0]
        num_reads = temp.total_mapped_reads_hmmcopy.values[0]
        
        temp_big = pd.DataFrame({
            'chr': [row.chr], 'start': [row.start], 'end': [row.end], 'reads': [reads], 'true_observed_CN': [true_observed_CN], 
            'true_rt_state': [rt_state], 'true_rt_norm': [rt_norm], 'A': [A], 'B': [B],
            'cell_id': [cell_id], 'gc_slope': [gc_slope], 'gc_int': [gc_int],
            'sigma1': [sigma1], 'gc': [row.gc], 'mcf7rt': [row.mcf7rt], 'true_G1_state': [state],
            'true_s_time': [true_s_time], 'true_frac_rt': [true_frac_rt], 'total_mapped_reads_hmmcopy': [num_reads]
        })
        
        aggregated_df.append(temp_big)
        
    aggregated_df = pd.concat(aggregated_df, ignore_index=True)
    aggregated_df['rpm'] = (aggregated_df['reads'] / sum(aggregated_df['reads'].values)) * 1e6
    
    return aggregated_df



def main():
    argv = get_args()

    ref_data = pd.read_csv(argv.ref_data, dtype={'Chromosome': str})
    ref_data = ref_data.dropna().reset_index(drop=True)
    ref_data.rename(columns={'Chromosome': 'chr', 'Start': 'start', 'End': 'end'}, inplace=True)
    ref_data['chr'] = ref_data['chr'].astype('category')
    
    ref_data_500kb = ref_data.loc[ref_data['bin_size']==500000]
    ref_data_500kb.reset_index(inplace=True, drop=True)
    print('ref_data_500kb.head()\n', ref_data_500kb.head())

    ref_data = ref_data.loc[ref_data['bin_size']==int(argv.bin_size)]
    ref_data.reset_index(inplace=True, drop=True)
    print('ref_data.head()\n', ref_data.head())

    gc_profile_500kb = pd.read_csv(argv.gc_profile_500kb)
    ref_data = filter_small_bin_df(gc_profile_500kb, ref_data)
    ref_data_500kb = filter_small_bin_df(gc_profile_500kb, ref_data_500kb)
    print('ref_data.head()\n', ref_data.head())
    print('ref_data_500kb.head()\n', ref_data_500kb.head())

    print('simulating S-phase cells...')
    diploid_s_cells = simulate_diploid_s_cells(
        ref_data, argv.num_cells_S, argv.gc_slope, argv.gc_int, argv.sigma1, argv.A, argv.B, argv.num_reads, 
        rt_col='mcf7rt', cell_id_prefix='cell_S', s_time_stdev=argv.s_time_stdev
    )
    print('diploid_s_cells.head()\n', diploid_s_cells.head())

    print('simulating G1-phase cells...')
    diploid_g_cells = simulate_diploid_g_cells(
        ref_data, argv.num_cells_G, argv.gc_slope, argv.gc_int, argv.sigma1, argv.A, argv.B, argv.num_reads, 
        rt_col='mcf7rt', cell_id_prefix='cell_G'
    )
    print('diploid_g_cells.head()\n', diploid_g_cells.head())

    # aggregate simulated data into 500kb bins if bin size was initially small
    # if argv.bin_size < 500000:
    #     print('aggregating S-phase cells into 500kb bins...')
    #     sim_s_to_500kb = []
    #     for cell_id, cell_cn in diploid_s_cells.groupby('cell_id'):
    #         aggregated_df = aggregate_small_bin_df(ref_data_500kb, cell_cn)
    #         sim_s_to_500kb.append(aggregated_df)
    #     diploid_s_cells = pd.concat(sim_s_to_500kb, ignore_index=True)
    #     print('diploid_s_cells.head()\n', diploid_s_cells.head())

    #     print('aggregating G1-phase cells into 500kb bins...')
    #     sim_g_to_500kb = []
    #     for cell_id, cell_cn in diploid_g_cells.groupby('cell_id'):
    #         aggregated_df = aggregate_small_bin_df(ref_data_500kb, cell_cn)
    #         sim_g_to_500kb.append(aggregated_df)
    #     diploid_g_cells = pd.concat(sim_g_to_500kb, ignore_index=True)
    #     print('diploid_g_cells.head()\n', diploid_g_cells.head())


    diploid_s_cells.to_csv(argv.s_out, sep='\t', index=False)
    diploid_g_cells.to_csv(argv.g_out, sep='\t', index=False)



if __name__ == '__main__':
    main()
