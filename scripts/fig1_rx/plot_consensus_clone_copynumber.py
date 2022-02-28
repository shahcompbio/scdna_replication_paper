import os
import matplotlib
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from scgenome.cnplot import plot_cell_cn_profile
from matplotlib.backends.backend_pdf import PdfPages
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('clone_states', help='input clone copynumber states tsv file')
    p.add_argument('clone_copy', help='input clone copynumber copy tsv file')
    p.add_argument('output_pdf', help='output pdf file showing consensus copy number profiles for each clone')

    return p.parse_args()


def plot_separate_pages(df, argv):
    ''' Plot each clone as a separate page with both copy and copy2_norm profiles. '''
    with PdfPages(argv.output_pdf) as pdf:
        for clone_id, plot_data in df.groupby('clone_id'):
            fig, ax = plt.subplots(1, 1, figsize=(16, 6))
            ax = ax.flatten()

            # plot copy colored by state
            _ = plot_cell_cn_profile(
                ax[0],
                plot_data,
                'copy',
                'state',
            )

            ax[0].set_title('clone {}'.format(clone_id))
            pdf.savefig(fig)
            plt.close()
            print('done plotting for clone {}'.format(clone_id))


def plot_as_one_page(df, argv, num_clones):
    ''' Plot only copy2_norm and don't use separate pages '''
    fig, ax = plt.subplots(num_clones, 1, figsize=(16, 6*num_clones))
    ax = ax.flatten()

    i = 0
    for clone_id, plot_data in df.groupby('clone_id'):
        _ = plot_cell_cn_profile(
            ax[i],
            plot_data,
            'copy',
            'state',
        )
        ax[i].set_title('clone {}'.format(clone_id))
        i += 1

    fig.savefig(argv.output_pdf, bbox_inches='tight')


def main():
    argv = get_args()
    clone_idx = ['chr', 'start', 'end']
    clone_states = pd.read_csv(argv.clone_states, sep='\t', index_col=clone_idx)
    clone_copy = pd.read_csv(argv.clone_copy, sep='\t', index_col=clone_idx)
    num_clones = len(clone_states.columns)

    # melt these three df into one df where the columns are 
    # chr, start, end, clone_id, state, copy, copy2_norm
    clone_states2 = pd.melt(clone_states.reset_index(),
                            id_vars=clone_idx, value_vars=clone_states.columns,
                            var_name='clone_id', value_name='state')
    clone_copy2 = pd.melt(clone_copy.reset_index(),
                          id_vars=clone_idx, value_vars=clone_copy.columns,
                          var_name='clone_id', value_name='copy')
    df = pd.merge(clone_states2, clone_copy2)

    # see if deleting dfs clears up memory to prevent segmentation fault
    del clone_states2
    del clone_copy2
    del clone_states
    del clone_copy

    # reset column types to be compatible with scgenome.cnplot
    df.loc[:, 'chr'] = df['chr'].astype('category')
    df.loc[:, 'state'] = df['state'].astype('category')
    df.loc[:, 'copy'] = df['copy'].astype('float')

    # plot_separate_pages(df, argv)

    plot_as_one_page(df, argv, num_clones)

    print('target output file: {}'.format(argv.output_pdf))
    print('target file exists: {}'.format(os.path.exists(argv.output_pdf)))
    return


if __name__ == '__main__':
    main()
    print('after main()')


