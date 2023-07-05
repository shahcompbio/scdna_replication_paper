import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import cross_validate
from sklearn.model_selection import RepeatedKFold
from argparse import ArgumentParser


def get_args():
    p = ArgumentParser()

    p.add_argument('embeddings', type=str, help='PCA embeddings with relevant per-profile features')
    p.add_argument('coefficients', type=str, help='table of regression coefficients for each clone feature')
    p.add_argument('plot', type=str, help='barplot of regression coefficients for each clone feature')

    return p.parse_args()



def main():
    argv = get_args()

    # load the clone embeddings
    clone_embeddings = pd.read_csv(argv.embeddings)

    # list of features to use in the regression
    # note that the signature and type columns are one-hot encoded
    feature_cols = [
        'ploidy', 'type_GM18507', 'type_T47D', 'type_hTERT', 'type_OV2295', 'type_HGSOC', 'type_TNBC',
        'signature_HRD', 'signature_TD', 'signature_FBI'
    ]

    cv = RepeatedKFold(n_splits=5, n_repeats=5, random_state=0)

    # transform the features to have mean 0 and variance 1
    X = StandardScaler().fit_transform(clone_embeddings[feature_cols])
    alphas = np.logspace(-10, 10, 21)
    feature_coef_df =[]
    for i in range(5):
        # get the ith PC
        ycol = 'PC{}'.format(i+1)
        y = clone_embeddings[ycol]
        # use the features to predict the PC on each fold
        cv_model = cross_validate(
            LinearRegression(),
            X,
            y,
            cv=cv,
            return_estimator=True,
            n_jobs=2,
        )
        # add the coefficients to the dataframe
        for est in cv_model["estimator"]:
            temp_df = pd.DataFrame({'feature': feature_cols, 'coef': est.coef_})
            temp_df['PC'] = ycol
            feature_coef_df.append(temp_df)

    feature_coef_df = pd.concat(feature_coef_df, ignore_index=True)

    # plot the log10 of absolute value of the coefficients
    feature_coef_df['log10_abs_coef'] = np.log10(abs(feature_coef_df['coef']))
    
    fig = plt.figure(figsize=(12, 4))
    sns.barplot(data=feature_coef_df, x='feature', y='log10_abs_coef', hue='PC')
    plt.ylabel('log10(abs(Coefficient))')
    plt.xlabel('Feature')
    plt.title('K-fold CV of linear regression coefficients for each PC')
    # reset the x-tick labels to rename '_' to '\n'
    plt.xticks(range(len(feature_cols)), [x.replace('_', '\n') for x in feature_cols])
    fig.savefig(argv.plot, dpi=300, bbox_inches='tight')

    # save the coefficients to a table
    feature_coef_df.to_csv(argv.coefficients, index=False)


if __name__ == '__main__':
    main()