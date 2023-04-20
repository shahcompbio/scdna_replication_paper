import scgenome
import pandas as pd

def main():
    genome_fasta = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/data/replication_timing/Homo_sapiens.GRCh37.70.dna.chromosomes.fa'
    prefix = '/juno/work/shah/users/weinera2/projects/scdna_replication_paper/data/replication_timing/'
    path_dict = {
        'mcf7_rt': 'wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig',
        'bg02es_rt': 'wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig',
        'bj_rt': 'wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig',
        'gm06990_rt': 'wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig',
        'gm12801_rt': 'wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig',
        'gm12812_rt': 'wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig',
        'gm12813_rt': 'wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig',
        'gm12878_rt': 'wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig',
        'helas3_rt': 'wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig',
        'hepg2_rt': 'wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig',
        'huvec_rt': 'wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig',
        'imr90_rt': 'wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig',
        'k562_rt': 'wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig',
        'sknsh_rt': 'wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig',
        'nhek_rt': 'wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig'
    }

    data = []
    # loop over different bin sizes
    for bin_size in (5e3, 1e4, 5e4, 5e5):
        bin_size = int(bin_size)
        print('bin size:', bin_size)

        # create bins at the desired size
        genome_binning = scgenome.tl.create_bins(bin_size)

        # count the gc content for the bins
        genome_binning = scgenome.tl.count_gc(genome_binning, genome_fasta, proportion=True)
        genome_df = genome_binning.as_df()

        # find all replication timing profiles at this bin size
        for cell_line, path in path_dict.items():
            print('cell line:', cell_line)
            temp_binning = scgenome.tl.mean_from_bigwig(genome_binning, prefix+path, cell_line, chr_prefix='chr')
            genome_df = pd.merge(genome_df, temp_binning.as_df())  # merge this cell line's RT as a new column
        
        # append bins of this size to list
        data.append(genome_df.assign(bin_size=bin_size))
    
    # concatenate all bin sizes into one df
    data = pd.concat(data, ignore_index=True)

    # save the gc and rt bin information to a csv file in the data folder
    data.to_csv('/juno/work/shah/users/weinera2/projects/scdna_replication_paper/data/gc_rt_bins.csv', index=False)


if __name__=='__main__':
    main()
