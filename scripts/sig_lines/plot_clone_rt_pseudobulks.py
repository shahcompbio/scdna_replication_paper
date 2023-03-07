import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from common.plot_utils import plot_cell_cn_profile2


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_s', type=str, help='full df for S-phase cells')
    p.add_argument('rt', type=str, help='table of RT pseudobulk profiles')
    p.add_argument('rep_col', type=str, help='column for replication state (relevant for RT pseudobulk profile column names)')
    p.add_argument('dataset', type=str, help='name of this dataset')
    p.add_argument('plot', type=str, help='figure comparing the clone pseudobulk rt profiles')

    return p.parse_args()


def plot_clone_rt_profiles(rt, argv):
    cols = [x for x in rt.columns if x.startswith('pseudobulk_clone') and x.endswith(argv.rep_col)]
    clones = [x.split('_')[1].replace('clone', '') for x in cols]
    ref_clone = clones[0]

    fig, ax = plt.subplots(4, 1, figsize=(16,16), tight_layout=True)
    ax = ax.flatten()
    i = 0
    for clone_id, col in zip(clones, cols):
        # plot the whole genome
        plot_cell_cn_profile2(
            ax[0], rt, col, color='C{}'.format(i), 
            max_cn=None, scale_data=False, lines=True, label=clone_id
        )
        # zoom in on chr1
        plot_cell_cn_profile2(
            ax[1], rt, col, color='C{}'.format(i), chromosome='1',
            max_cn=None, scale_data=False, lines=True, label=clone_id
        )
        
        # compute distance between this clone's RT and the reference clone (A)
        if clone_id != ref_clone:
            rt['temp_diff_rt'] = rt[col] - rt['pseudobulk_clone{}_{}'.format(ref_clone, argv.rep_col)]

            # plot the whole genome
            plot_cell_cn_profile2(
                ax[2], rt, 'temp_diff_rt', color='C{}'.format(i), 
                max_cn=None, scale_data=False, lines=True, label=clone_id
            )
            # zoom in on chr1
            plot_cell_cn_profile2(
                ax[3], rt, 'temp_diff_rt', color='C{}'.format(i), chromosome='1',
                max_cn=None, scale_data=False, lines=True, label=clone_id
            )
            
        i += 1

    for i in range(4):
        ax[i].set_title(argv.dataset)
        ax[i].legend(title='Clone ID')
        if i < 2:
            ax[i].set_ylabel('RT profile\n<--late | early-->')
        else:
            ax[i].set_ylabel('Relative RT to clone {}\n<--clone later | clone earlier-->'.format(ref_clone))

    fig.savefig(argv.plot, bbox_inches='tight', dpi=300)


def main():
    argv = get_args()
    # load long-form dataframe from different cell cycle phases
    cn_s = pd.read_csv(argv.cn_s, sep='\t')

    # load pseudobulk RT profiles
    rt = pd.read_csv(argv.rt, sep='\t')

    cn_s.chr = cn_s.chr.astype('str')
    cn_s.chr = cn_s.chr.astype('category')
    rt.chr = rt.chr.astype('str')
    rt.chr = rt.chr.astype('category')

    # merge end and gc columns into RT pseudobulks
    rt = pd.merge(rt, cn_s[['chr', 'start', 'end', 'gc']].drop_duplicates())

    # plot pseudoublk RT profiles
    plot_clone_rt_profiles(rt, argv)


if __name__ == '__main__':
    main()
