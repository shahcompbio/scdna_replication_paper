import pandas as pd
import numpy as np
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('-i', '--input', type=str, nargs='+', help='cell cycle clone counts for each dataset')
    p.add_argument('-o', '--output', type=str, help='summary of clone counts for each dataset in one tsv file')

    return p.parse_args()


def load_clone_counts(argv):
    """ Use the  cell cycle clone counts input to create a table of the number of cells in each clone for each cell cycle phase and dataset. """
    counts = []
    for path in argv.input:
        # load the counts
        temp_counts = pd.read_csv(path, sep='\t')
        # subset to just the columns with clone-id and cell counts
        temp_counts = temp_counts[['clone_id', 'num_cells_s', 'num_cells_g']]
        # sum num_cells_s column across all rows with the same clone_id
        # this is important as there are multiple rows (libraries) for each clone_id
        temp_counts = temp_counts.groupby('clone_id').sum().reset_index()
        # add dataset name as column
        d = path.split('/')[2]
        temp_counts['dataset'] = d
        # compute the fraction of cells in each clone for each cell cycle phase
        temp_counts['frac_cells_s'] = temp_counts['num_cells_s'] / (temp_counts['num_cells_s'] + temp_counts['num_cells_g'])
        temp_counts['frac_cells_g'] = temp_counts['num_cells_g'] / (temp_counts['num_cells_s'] + temp_counts['num_cells_g'])
        # add to list of counts
        counts.append(temp_counts)
    # concatenate all counts
    counts = pd.concat(counts, ignore_index=True)

    return counts


def main():
    argv = get_args()
    counts = load_clone_counts(argv)
    counts.to_csv(argv.output, sep='\t', index=False)


if __name__ == '__main__':
    main()
