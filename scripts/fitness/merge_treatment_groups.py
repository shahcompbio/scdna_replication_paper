import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('SA1035U_in', type=str, help='CN data from SA1035 untreated sample')
    p.add_argument('SA1035T_in', type=str, help='CN data from SA1035 treated sample')
    p.add_argument('SA535U_in', type=str, help='CN data from SA535 untreated sample')
    p.add_argument('SA535T_in', type=str, help='CN data from SA535 treated sample')
    p.add_argument('SA1035_out', type=str, help='merged CN data from all SA1035 cells')
    p.add_argument('SA535_out', type=str, help='merged CN data from all SA535 cells')

    return p.parse_args()


def main():
    argv = get_args()

    # load in data from both cell cycle phases and treatment groups
    cn_SA1035U = pd.read_csv(argv.SA1035U_in, sep='\t')
    cn_SA1035T = pd.read_csv(argv.SA1035T_in, sep='\t')
    cn_SA535U = pd.read_csv(argv.SA535U_in, sep='\t')
    cn_SA535T = pd.read_csv(argv.SA535T_in, sep='\t')

    # label which cells are treated vs untreated so I can split later
    cn_SA1035U['treated'] = False
    cn_SA535U['treated'] = False
    cn_SA1035T['treated'] = True
    cn_SA535T['treated'] = True

    # combine both treated and untreated cells into one dataframe
    cn_SA1035 = pd.concat([cn_SA1035T, cn_SA1035U], ignore_index=True)
    cn_SA535 = pd.concat([cn_SA535T, cn_SA535U], ignore_index=True)

    # save the merged samples for each PDX
    cn_SA1035.to_csv(argv.SA1035_out, sep='\t', index=False)
    cn_SA535.to_csv(argv.SA535_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
