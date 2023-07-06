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

    target_cols = [
        'PC1', 'PC2', 'PC3', 'PC4', 'PC5'
    ]

    cv = RepeatedKFold(n_splits=5, n_repeats=5, random_state=0)

    # transform the features to have mean 0 and variance 1
    X = StandardScaler().fit_transform(clone_embeddings[feature_cols])
    y = clone_embeddings[target_cols]

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
    feature_coef_df =[]
    for est in cv_model["estimator"]:
        for i, ycol in enumerate(target_cols):
            temp_df = pd.DataFrame({'feature': feature_cols, 'coef': est.coef_[i]})
            temp_df['PC'] = ycol
            feature_coef_df.append(temp_df)

    feature_coef_df = pd.concat(feature_coef_df, ignore_index=True)

    # plot the log10 of absolute value of the coefficients
    feature_coef_df['log10_abs_coef'] = np.log10(abs(feature_coef_df['coef']))

    # remove rows from feature_coef_df where coef==0 as these represent cases where the cross-validation resulted in a model that did not use that feature
    feature_coef_df = feature_coef_df.query('coef!=0')

    # create a new column named 'metafeature' which merges all the one-hot encoded features into a single feature
    feature_coef_df['metafeature'] = feature_coef_df['feature'].apply(lambda x: x.split('_')[0])
    
    fig, ax = plt.subplots(1, 4, figsize=(16, 4))
    # merge the first three subplots into one supblot
    ax[0].remove()
    ax[1].remove()
    ax[2].remove()
    # create a new subplot in the first position that spans the first three subplots
    ax[0] = fig.add_subplot(1, 4, (1, 3))

    sns.barplot(data=feature_coef_df, x='feature', y='log10_abs_coef', hue='PC', ax=ax[0])
    ax[0].set_ylabel('log10(abs(Coefficient))')
    ax[0].set_xlabel('Feature')
    ax[0].set_title('Multivariate linear regression of PCs\nusing K-fold CV')
    # reset the x-tick labels to rename '_' to '\n'
    ax[0].set_xticklabels([x.replace('_', '\n') for x in feature_cols])

    sns.barplot(data=feature_coef_df, x='metafeature', y='log10_abs_coef', ax=ax[3])
    ax[3].set_ylabel('log10(abs(Coefficient))')
    ax[3].set_xlabel('Feature')
    ax[3].set_title('Multivariate linear regression of PCs\nusing K-fold CV')

    # save the figure
    fig.savefig(argv.plot, bbox_inches='tight', dpi=300)

    # save the coefficients to a table
    feature_coef_df.to_csv(argv.coefficients, index=False)


if __name__ == '__main__':
    main()