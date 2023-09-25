import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('input', type=str, help='table of loss values per iteration')
    p.add_argument('output', type=str, help='output plot of loss curve')

    return p.parse_args()


def main():
    argv = get_args()

    df = pd.read_csv(argv.input)

    fig = plt.figure(figsize=(5, 2))
    plt.plot(df['iteration'], df['loss'])
    plt.xlabel("SVI step")
    plt.ylabel("ELBO loss")
    plt.title('Modeling clone RT profiles')
    
    fig.savefig(argv.output, bbox_inches='tight', dpi=300)


if __name__ == '__main__':
    main()
