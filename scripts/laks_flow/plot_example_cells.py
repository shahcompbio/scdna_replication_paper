import pandas as pd
import matplotlib.pyplot as plt
# save svg text as text, not curves
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['pdf.use14corefonts'] = True
import seaborn as sns
from scdna_replication_tools.plot_utils import plot_cell_cn_profile2
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('cn_t_g1_in', help='long-form copy number dataframe for T47D cells in G1 phase')
    p.add_argument('cn_t_g2_in', help='long-form copy number dataframe for T47D cells in G2 phase')
    p.add_argument('cn_t_s_in', help='long-form copy number dataframe for T47D cells in S phase')
    p.add_argument('cn_gm_g1_in', help='long-form copy number dataframe for GM18507 cells in G1 phase')
    p.add_argument('cn_gm_g2_in', help='long-form copy number dataframe for GM18507 cells in G2 phase')
    p.add_argument('cn_gm_s_in', help='long-form copy number dataframe for GM18507 cells in S phase')
    p.add_argument('cn_t_g1_out', help='representative cell profiles for T47D cells in G1 phase')
    p.add_argument('cn_t_g2_out', help='representative cell profiles for T47D cells in G2 phase')
    p.add_argument('cn_t_s_out', help='representative cell profiles for T47D cells in S phase')
    p.add_argument('cn_gm_g1_out', help='representative cell profiles for GM18507 cells in G1 phase')
    p.add_argument('cn_gm_g2_out', help='representative cell profiles for GM18507 cells in G2 phase')
    p.add_argument('cn_gm_s_out', help='representative cell profiles for GM18507 cells in S phase')

    return p.parse_args()


def plot_cell_cn_profile_with_density(cell_cn, title=None):
    # create figure of shape 15 x 3
    fig = plt.figure(figsize=(8, 2))

    # add subplot for cell cn profile
    # this plot should have a height of 1, width of 0.9, and be placed at (0.0, 0.0)
    ax = fig.add_axes([0.0, 0.0, 0.9, 1.0])
    plot_cell_cn_profile2(ax, cell_cn, 'rpm', cn_field_name='state', rawy=True)
    ax.set_ylabel('reads per million', fontname='Arial', fontsize=11)
    # set the x-axis label to be Arial of size 11
    ax.set_xlabel('chromosome', fontname='Arial', fontsize=11)

    if title is not None:
        ax.set_title(title, fontname='Arial', fontsize=11)
    
    # Set the font name for axis tick labels to be Arial
    for tick in ax.get_xminorticklabels():
        tick.set_fontname('Arial')
        tick.set_fontsize(8)
    for tick in ax.get_yticklabels():
        tick.set_fontname('Arial')
        tick.set_fontsize(8)


    # add subplot for density plot
    # this plot should have a height of 1, width of 0.1, and be placed at (0.9, 0.0)
    ax = fig.add_axes([0.9, 0.0, 0.1, 1.0])
    sns.kdeplot(data=cell_cn, y='rpm', ax=ax, shade=True, color='k')
    ax.set_yticks([])
    ax.set_xticks([])
    ax.set_ylabel('')
    ax.set_xlabel('')
    # turn off the spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    return fig


def main():
    argv = get_args()

    # load in all the input dataframes
    cn_t_g1 = pd.read_csv(argv.cn_t_g1_in, sep='\t')
    cn_t_g2 = pd.read_csv(argv.cn_t_g2_in, sep='\t')
    cn_t_s = pd.read_csv(argv.cn_t_s_in, sep='\t')
    cn_gm_g1 = pd.read_csv(argv.cn_gm_g1_in, sep='\t')
    cn_gm_g2 = pd.read_csv(argv.cn_gm_g2_in, sep='\t')
    cn_gm_s = pd.read_csv(argv.cn_gm_s_in, sep='\t')

    # convert chr to category type
    cn_t_g1['chr'] = cn_t_g1['chr'].astype('category')
    cn_t_g2['chr'] = cn_t_g2['chr'].astype('category')
    cn_t_s['chr'] = cn_t_s['chr'].astype('category')
    cn_gm_g1['chr'] = cn_gm_g1['chr'].astype('category')
    cn_gm_g2['chr'] = cn_gm_g2['chr'].astype('category')
    cn_gm_s['chr'] = cn_gm_s['chr'].astype('category')


    # plot a T47D cell in G1 phase
    for cell_id, cell_cn in cn_t_g1.groupby('cell_id'):
        fig = plot_cell_cn_profile_with_density(cell_cn, title='{}\nT47D, FACS G1-phase'.format(cell_id))
        fig.savefig(argv.cn_t_g1_out, dpi=300, bbox_inches='tight')
        break

    # plot a T47D cell in G2 phase
    for cell_id, cell_cn in cn_t_g2.groupby('cell_id'):
        fig = plot_cell_cn_profile_with_density(cell_cn, title='{}\nT47D, FACS G2-phase'.format(cell_id))
        fig.savefig(argv.cn_t_g2_out, dpi=300, bbox_inches='tight')
        break

    # plot a T47D cell in S phase
    for cell_id, cell_cn in cn_t_s.groupby('cell_id'):
        fig = plot_cell_cn_profile_with_density(cell_cn, title='{}\nT47D, FACS S-phase'.format(cell_id))
        fig.savefig(argv.cn_t_s_out, dpi=300, bbox_inches='tight')
        break

    # plot a GM18507 cell in G1 phase
    for cell_id, cell_cn in cn_gm_g1.groupby('cell_id'):
        fig = plot_cell_cn_profile_with_density(cell_cn, title='{}\nGM18507, FACS G1-phase'.format(cell_id))
        fig.savefig(argv.cn_gm_g1_out, dpi=300, bbox_inches='tight')
        break

    # plot a GM18507 cell in G2 phase
    for cell_id, cell_cn in cn_gm_g2.groupby('cell_id'):
        fig =  plot_cell_cn_profile_with_density(cell_cn, title='{}\nGM18507, FACS G2-phase'.format(cell_id))
        fig.savefig(argv.cn_gm_g2_out, dpi=300, bbox_inches='tight')
        break

    # plot a GM18507 cell in S phase
    for cell_id, cell_cn in cn_gm_s.groupby('cell_id'):
        fig = plot_cell_cn_profile_with_density(cell_cn, title='{}\nGM18507, FACS S-phase'.format(cell_id))
        fig.savefig(argv.cn_gm_s_out, dpi=300, bbox_inches='tight')
        break



if __name__ == '__main__':
    main()
