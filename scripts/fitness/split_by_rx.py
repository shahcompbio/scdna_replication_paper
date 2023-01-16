import pandas as pd
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('SA1035_s_in', type=str, help='S-phase CN data from SA1035 treated and untreated samples')
    p.add_argument('SA1035_g_in', type=str, help='G1/2-phase CN data from SA1035 treated and untreated samples')
    p.add_argument('SA535_s_in', type=str, help='S-phase CN data from SA535 treated and untreated samples')
    p.add_argument('SA535_g_in', type=str, help='G1/2-phase CN data from SA535 treated and untreated samples')
    p.add_argument('SA609_s_in', type=str, help='S-phase CN data from SA609 treated and untreated samples')
    p.add_argument('SA609_g_in', type=str, help='G1/2-phase CN data from SA609 treated and untreated samples')
    p.add_argument('SA1035U_s_out', type=str, help='S-phase CN data from SA1035 untreated sample')
    p.add_argument('SA1035U_g_out', type=str, help='G1/2-phase CN data from SA1035 untreated sample')
    p.add_argument('SA1035T_s_out', type=str, help='S-phase CN data from SA1035 treated sample')
    p.add_argument('SA1035T_g_out', type=str, help='G1/2-phase CN data from SA1035 treated sample')
    p.add_argument('SA535U_s_out', type=str, help='S-phase CN data from SA535 untreated sample')
    p.add_argument('SA535U_g_out', type=str, help='G1/2-phase CN data from SA535 untreated sample')
    p.add_argument('SA535T_s_out', type=str, help='S-phase CN data from SA535 treated sample')
    p.add_argument('SA535T_g_out', type=str, help='G1/2-phase CN data from SA535 treated sample')
    p.add_argument('SA609U_s_out', type=str, help='S-phase CN data from SA609 untreated sample (line 1)')
    p.add_argument('SA609U_g_out', type=str, help='G1/2-phase CN data from SA609 untreated sample (line 1)')
    p.add_argument('SA609T_s_out', type=str, help='S-phase CN data from SA609 treated sample (line 1)')
    p.add_argument('SA609T_g_out', type=str, help='G1/2-phase CN data from SA609 treated sample (line 1)')
    p.add_argument('SA609U2_s_out', type=str, help='S-phase CN data from SA609 untreated sample (line 2)')
    p.add_argument('SA609U2_g_out', type=str, help='G1/2-phase CN data from SA609 untreated sample (line 2)')
    p.add_argument('SA609T2_s_out', type=str, help='S-phase CN data from SA609 treated sample (line 2)')
    p.add_argument('SA609T2_g_out', type=str, help='G1/2-phase CN data from SA609 treated sample (line 2)')
    

    return p.parse_args()


def main():
    argv = get_args()

    # load in data from both cell cycle phases
    cn_s_SA1035 = pd.read_csv(argv.SA1035_s_in, sep='\t')
    cn_g_SA1035 = pd.read_csv(argv.SA1035_g_in, sep='\t')
    cn_s_SA535 = pd.read_csv(argv.SA535_s_in, sep='\t')
    cn_g_SA535 = pd.read_csv(argv.SA535_g_in, sep='\t')
    cn_s_SA609 = pd.read_csv(argv.SA609_s_in, sep='\t')
    cn_g_SA609 = pd.read_csv(argv.SA609_g_in, sep='\t')

    # split each dataframe into treated and untreated subsets
    cn_s_SA1035U = cn_s_SA1035.loc[cn_s_SA1035['treated']==False]
    cn_g_SA1035U = cn_g_SA1035.loc[cn_g_SA1035['treated']==False]
    cn_s_SA1035T = cn_s_SA1035.loc[cn_s_SA1035['treated']==True]
    cn_g_SA1035T = cn_g_SA1035.loc[cn_g_SA1035['treated']==True]
    cn_s_SA535U = cn_s_SA535.loc[cn_s_SA535['treated']==False]
    cn_g_SA535U = cn_g_SA535.loc[cn_g_SA535['treated']==False]
    cn_s_SA535T = cn_s_SA535.loc[cn_s_SA535['treated']==True]
    cn_g_SA535T = cn_g_SA535.loc[cn_g_SA535['treated']==True]

    # split SA609 dataframes by both treatment and line
    cn_s_SA609U = cn_s_SA609.loc[(cn_s_SA609['treated']==False) & (cn_s_SA609['line']==1)]
    cn_g_SA609U = cn_g_SA609.loc[(cn_g_SA609['treated']==False) & (cn_g_SA609['line']==1)]
    cn_s_SA609T = cn_s_SA609.loc[(cn_s_SA609['treated']==True) & (cn_s_SA609['line']==1)]
    cn_g_SA609T = cn_g_SA609.loc[(cn_g_SA609['treated']==True) & (cn_g_SA609['line']==1)]
    cn_s_SA609U2 = cn_s_SA609.loc[(cn_s_SA609['treated']==False) & (cn_s_SA609['line']==2)]
    cn_g_SA609U2 = cn_g_SA609.loc[(cn_g_SA609['treated']==False) & (cn_g_SA609['line']==2)]
    cn_s_SA609T2 = cn_s_SA609.loc[(cn_s_SA609['treated']==True) & (cn_s_SA609['line']==2)]
    cn_g_SA609T2 = cn_g_SA609.loc[(cn_g_SA609['treated']==True) & (cn_g_SA609['line']==2)]

    # save samples that are now split by both cell cycle and treatement status
    cn_s_SA1035U.to_csv(argv.SA1035U_s_out, sep='\t', index=False)
    cn_g_SA1035U.to_csv(argv.SA1035U_g_out, sep='\t', index=False)
    cn_s_SA1035T.to_csv(argv.SA1035T_s_out, sep='\t', index=False)
    cn_g_SA1035T.to_csv(argv.SA1035T_g_out, sep='\t', index=False)
    cn_s_SA535U.to_csv(argv.SA535U_s_out, sep='\t', index=False)
    cn_g_SA535U.to_csv(argv.SA535U_g_out, sep='\t', index=False)
    cn_s_SA535T.to_csv(argv.SA535T_s_out, sep='\t', index=False)
    cn_g_SA535T.to_csv(argv.SA535T_g_out, sep='\t', index=False)
    cn_s_SA609U.to_csv(argv.SA609U_s_out, sep='\t', index=False)
    cn_g_SA609U.to_csv(argv.SA609U_g_out, sep='\t', index=False)
    cn_s_SA609T.to_csv(argv.SA609T_s_out, sep='\t', index=False)
    cn_g_SA609T.to_csv(argv.SA609T_g_out, sep='\t', index=False)
    cn_s_SA609U2.to_csv(argv.SA609U2_s_out, sep='\t', index=False)
    cn_g_SA609U2.to_csv(argv.SA609U2_g_out, sep='\t', index=False)
    cn_s_SA609T2.to_csv(argv.SA609T2_s_out, sep='\t', index=False)
    cn_g_SA609T2.to_csv(argv.SA609T2_g_out, sep='\t', index=False)


if __name__ == '__main__':
    main()
