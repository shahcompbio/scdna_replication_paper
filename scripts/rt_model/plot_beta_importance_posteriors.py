import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='table of beta importance posteriors')
    p.add_argument('--noX', action='store_true', help='binary indicator to include for version of the model where chrX is excluded')
    p.add_argument('output', type=str, help='histogram of all the beta importance posteriors')

    return p.parse_args()


def main():
    argv = get_args()

    # read in the table of beta importance posteriors
    beta_importance = pd.read_csv(argv.input)

    # plot the histogram of beta importance posteriors
    # each column has its own histogram and color
    # the x-axis is the relative importance of each beta
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, col in enumerate(beta_importance.columns):
        # remove the beta_ prefix and _importance suffix from the column name for plotting
        label = col.replace('beta_', '').replace('_importance', '').replace('_', ' ')
        sns.distplot(beta_importance[col], hist=False, kde=True, color=sns.color_palette()[i], kde_kws={'linewidth': 2}, label=col)
    ax.set_xlabel('Beta coefficient\n<-- less important | more important -->')
    ax.set_ylabel('Density')

    if argv.noX:
        ax.set_title('Feature importance\n(no chrX)')
    else:
        ax.set_title('Feature importance')
    
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)

    fig.savefig(argv.output, dpi=300, bbox_inches='tight')


if __name__ == '__main__':
    main()
