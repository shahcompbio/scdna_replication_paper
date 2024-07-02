import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='table of scRNA cell cycle counts for each site')
    p.add_argument('plot', type=str, help='figure showing the distribution of cell cycle phases for each site')

    return p.parse_args()


def main():
    argv = get_args()

    # load the table
    df = pd.read_csv(argv.input)

    # count the number of cells within each site & phase combination
    count_df = df[['site', 'Phase']].value_counts().reset_index().rename(columns={0: 'num_cells', 'Phase': 'phase'}).sort_values(by=['site', 'phase']).reset_index(drop=True)

    # plot the distribution of cell cycle phases for each site as pie charts
    fig, ax = plt.subplots(1, 2, figsize=(8, 4), tight_layout=True)

    ax[0].pie(count_df.query('site=="LEFT_ADNEXA"')['num_cells'], labels=count_df.query('site=="LEFT_ADNEXA"')['phase'], autopct='%1.1f%%', shadow=True, startangle=90)
    ax[1].pie(count_df.query('site=="INFRACOLIC_OMENTUM"')['num_cells'], labels=count_df.query('site=="INFRACOLIC_OMENTUM"')['phase'], autopct='%1.1f%%', shadow=True, startangle=90)

    ax[0].set_title('Left Adnexa (primary)\nNGD clones dominant')
    ax[1].set_title('Infracolic Omentum (met)\nWGD clone dominant')
    fig.suptitle('scRNA cell cycle distribution of tumor cells')

    # save the figure
    fig.savefig(argv.plot, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
