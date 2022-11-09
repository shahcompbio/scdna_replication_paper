from argparse import ArgumentParser
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation


def get_args():
    p = ArgumentParser()

    p.add_argument('-ir', '--input_rt', type=str, nargs='+', help='pseduobulk rt profiles for each dataset')
    p.add_argument('-ic', '--input_cn', type=str, nargs='+', help='pseduobulk cn profiles for each dataset')
    p.add_argument('-r', '--rep_col', type=str, help='column used to denote replication state')
    p.add_argument('-d', '--dataset', type=str, nargs='+', help='list of dataset names')
    p.add_argument('--table', help='table of all the concatenated input tables')
    p.add_argument('--plot', help='Summary plots of subclonal RT diffs across all datasets')

    return p.parse_args()


def load_rt_data(argv):
    # load subclonal rt diff table for each dataset
    rt = []
    for rt_path, d in zip(argv.input_rt, argv.dataset):
        temp_rt = pd.read_csv(rt_path, sep='\t')
        rt.append(temp_rt)

    # concatenate into one df
    rt = pd.concat(rt, ignore_index=True)
    return rt