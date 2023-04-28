import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from argparse import ArgumentParser

def get_args():
    p = ArgumentParser()

    p.add_argument('input', help='table of S and G1/2-phase cell counts for each clone in the whole cohort')
    p.add_argument('output', help='figure showing S-phase fraction for each clone in the whole cohort')

    return p.parse_args()


# TODO: verify that this bootstrap method is correct
def compute_spf_with_bootstrapping(df, n_boot=100, boot_frac=0.75):
    '''
    For each row in df, compute the S-phase fraction (spf) and the 95% confidence interval of the spf using bootstrapping.
    '''
    spf = []
    spf_low = []
    spf_high = []
    for i in range(len(df)):
        # get the number of cells in S-phase and G1/2-phase for this clone
        num_cells_s = df.iloc[i]['num_cells_s']
        num_cells_g = df.iloc[i]['num_cells_g']
        total_cells = num_cells_s + num_cells_g
        temp_spf = num_cells_s / total_cells
        # compute the spf for this clone
        spf.append(temp_spf)
        # bootstrap the number of cells in S-phase and G1/2-phase for this clone using just a fraction of the data
        spf_boot = []
        for j in range(n_boot):
            # sample a bernoulli random variable with probability of success equal to the spf
            bernoulli_flips = np.random.binomial(total_cells, temp_spf)
            # compute the spf for this bootstrapped sample
            spf_boot.append(bernoulli_flips / total_cells)

        # compute the 95% confidence interval of the spf
        spf_low.append(np.percentile(spf_boot, 2.5))
        spf_high.append(np.percentile(spf_boot, 97.5))
    # convert the lists to numpy arrays
    spf = np.array(spf)
    spf_yerr = np.array([spf - spf_low, spf_high - spf])
    return spf, spf_yerr


def main():
    argv = get_args()

    df = pd.read_csv(argv.input)

   # compute the fraction of cells in each clone that are in S-phase 
    # df['spf'] = df['num_cells_s'] / (df['num_cells_s'] + df['num_cells_g'])
    df['spf'], spf_yerr = compute_spf_with_bootstrapping(df, n_boot=100)

    # plot a barplot of the S-phase fraction for each clone with a legend on the x-axis that denotes the dataset
    dataset_to_color = {}
    cmap = plt.get_cmap('tab20')
    for i, dataset in enumerate(df['dataset'].unique()):
        dataset_to_color[dataset] = cmap(i/len(df['dataset'].unique()))
    df['color'] = df['dataset'].apply(lambda x: dataset_to_color[x])

    # the ID for each row is the dataset + clone ID
    df['id'] = df['dataset'] + ' clone ' + df['clone_id']

    fig, ax = plt.subplots(2, 1, figsize=(12, 8), tight_layout=True)
    
    # plot the S-phase fraction for each clone
    for dataset in df['dataset'].unique():
        idx = df[df['dataset'] == dataset].index
        ax[0].bar(x=idx, height=df[df['dataset'] == dataset]['spf'], color=dataset_to_color[dataset], label=dataset, yerr=spf_yerr[:, idx])
    ax[0].set_xticklabels('')
    ax[0].set_ylabel('S-phase fraction (95% CI)')
    ax[0].set_title('All clones in the cohort')
    ax[0].set_xlim(-1, len(df))
    ax[0].legend(loc='upper right', ncol=2)

    # plot the S-phase fraction for each sample
    for i, dataset in enumerate(df['dataset'].unique()):
        temp_num_s = df[df['dataset'] == dataset]['num_cells_s'].sum()
        temp_num_g = df[df['dataset'] == dataset]['num_cells_g'].sum()
        temp_total_cells = temp_num_s + temp_num_g
        temp_spf = temp_num_s / temp_total_cells
        # do a boostrap to get the 95% confidence interval of the spf
        spf_boot = []
        for j in range(100):
            bernoulli_flips = np.random.binomial(temp_total_cells, temp_spf)
            spf_boot.append(bernoulli_flips / temp_total_cells)
        spf_low = np.percentile(spf_boot, 2.5)
        spf_high = np.percentile(spf_boot, 97.5)
        temp_yerr = np.array([[temp_spf - spf_low], [spf_high - temp_spf]])
        print('temp_yerr', temp_yerr)
        print('temp_spf', temp_spf)
        ax[1].bar(x=i, height=temp_spf, color=dataset_to_color[dataset], label=dataset, yerr=temp_yerr)
    ax[1].set_xticklabels('')
    ax[1].set_ylabel('S-phase fraction (95% CI)')
    ax[1].set_title('All samples in the cohort')
    ax[1].set_xlim(-1, len(df['dataset'].unique()))
    ax[1].legend(loc='upper right', ncol=2)
        
    fig.savefig(argv.output, dpi=300)



if __name__ == '__main__':
    main()
