'''
Created on Feb 21, 2023
@author: adamcweiner
'''
# from __future__ import division

import argparse
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats.mstats import mquantiles
from statsmodels.nonparametric.smoothers_lowess import lowess


class CorrectReadCount(object):
    """
    fit lowess/polynomial curve smoothing to reads-gc
    use fitted model to predict corrected gc and mappability
    values
    """

    def __init__(self, input, output, mappability=0.9,
                 smoothing_function='lowess',
                 polynomial_degree=2):
        self.mappability = mappability

        self.input = input
        self.output = output

    def valid(self, df):
        """adds valid column (calls with atleast one reads and non negative gc)
        :params df: pandas dataframe
        """

        df.loc[:, "valid"] = True

        df.loc[(df["reads"] <= 0) | (df['gc'] < 0), "valid"] = False

        return df

    def ideal(self, df):
        """adds ideal column
        :params df: pandas dataframe
        """
        df.loc[:, "ideal"] = True

        valid_reads = df[df["valid"]]["reads"]
        valid_gc = df[df["valid"]]["gc"]

        routlier = 0.01
        doutlier = 0.001

        range_l, range_h = mquantiles(valid_reads, prob=[0, 1 - routlier],
                                      alphap=1, betap=1)
        domain_l, domain_h = mquantiles(valid_gc, prob=[doutlier, 1 - doutlier],
                                        alphap=1, betap=1)

        df.loc[(df["valid"] == False) |
               (df["map"] < self.mappability) |
               (df["reads"] <= range_l) |
               (df["reads"] > range_h) |
               (df["gc"] < domain_l) |
               (df["gc"] > domain_h),
               "ideal"] = False

        return df

    def modal_quantile_regression(self, df_regression, lowess_frac=0.2):
        '''
        Compute quantile regression curves and select the modal quantile.
        '''
        # 2nd order polynomial quantile regression 

        q_range = range(10, 91, 1)
        quantiles = np.array(q_range) / 100
        quantile_names = [str(x) for x in q_range]

        # need at least 3 values to compute the quantiles
        if len(df_regression) < 10:
            return df_regression

        poly2_quantile_model = smf.quantreg('reads ~ gc + I(gc ** 2.0)', data=df_regression)
        poly2_quantile_fit = [poly2_quantile_model.fit(q=q) for q in quantiles]
        poly2_quantile_predict = [poly2_quantile_fit[i].predict(df_regression) for i in range(len(quantiles))]

        poly2_quantile_params = pd.DataFrame()

        for i in range(len(quantiles)):
            df_regression[quantile_names[i]] = poly2_quantile_predict[i]
            poly2_quantile_params[quantile_names[i]] = poly2_quantile_fit[i].params

        # integration and mode selection

        gc_min = df_regression['gc'].quantile(q=0.10)
        gc_max = df_regression['gc'].quantile(q=0.90)

        poly2_quantile_integration = np.zeros(len(quantiles) + 1)

        for i in range(len(quantiles)):
            params = poly2_quantile_params[quantile_names[i]].tolist()
            params.reverse()
            poly2 = np.poly1d(params)
            integ = poly2.integ()
            integrand = integ(gc_max) - integ(gc_min)
            poly2_quantile_integration[i + 1] = integrand

        # find the modal quantile

        distances = poly2_quantile_integration[1:] - poly2_quantile_integration[:-1]

        df_dist = pd.DataFrame({'quantiles': quantiles, 'quantile_names': quantile_names, 'distances': distances})
        dist_max = df_dist['distances'].quantile(q=0.95)
        df_dist_filter = df_dist[df_dist['distances'] < dist_max]
        df_dist_filter['lowess'] = lowess(df_dist_filter['distances'], df_dist_filter['quantiles'], frac=lowess_frac,
                                          return_sorted=False)

        modal_quantile = df_dist_filter.set_index('quantile_names')['lowess'].idxmin()

        # add values to table

        df_regression['modal_quantile'] = modal_quantile
        df_regression['modal_curve'] = df_regression[modal_quantile]
        df_regression['modal_corrected'] = df_regression['reads'] / df_regression[modal_quantile]

        return df_regression

    def write(self, df):
        """write results to the output file
        :param df: pandas dataframe
        """

        df.to_csv(self.output, index=False, sep=',', na_rep="NA")

    def main(self):
        # read input dataframe containing columns 
        # 'chr', 'start', 'end', 'width', 'gc', 'map', 'reads', 'cell_id'
        input_df = pd.read_csv(self.input)

        # store each cell's dataframe in a list of dataframes
        output_df = []

        for cell_id, df in input_df.groupby('cell_id'):

            # annotate valid and ideal bins
            df = self.valid(df)
            df = self.ideal(df)

            df['modal_quantile'] = 'NaN'
            df['modal_curve'] = 'NaN'
            df['modal_corrected'] = 'NaN'

            # filtering and sorting
            df_valid_gc = df[df['gc'] > 0]

            df_non_zero = df_valid_gc[df_valid_gc['reads'] > 0]

            df_regression = pd.DataFrame.copy(df_non_zero)

            df_regression.sort_values(by='gc', inplace=True)

            # modal quantile regression
            df_regression = self.modal_quantile_regression(df_regression, lowess_frac=0.2)

            # map results back to full data frame
            df.ix[df_regression.index, 'modal_quantile'] = df_regression['modal_quantile']
            df.ix[df_regression.index, 'modal_curve'] = df_regression['modal_curve']
            df.ix[df_regression.index, 'modal_corrected'] = df_regression['modal_corrected']

            # filter by mappability
            df['copy'] = df['modal_corrected']
            df['copy'][df['map'] < self.mappability] = float('NaN')

            df = df.rename(columns=({"modal_corrected": "cor_gc"}))

            df["cor_map"] = float("NaN")

            output_df.append(df)
    
        # concatenate the per-cell dataframes into one output dataframe
        output_df = pd.concat(output_df, ignore_index=True)

        # save
        self.write(output_df)


def parse_args():
    """
    parses command line arguments
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('input',
                        help='path to the input read count csv'
                        )

    parser.add_argument('output',
                        help='path to the output csv file'
                        )

    parser.add_argument('--mappability',
                        default=0.9,
                        type=float,
                        help='specify mappability threshold')

    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = parse_args()

    corr = CorrectReadCount(args.input, args.output,
                            mappability=args.mappability,
                            )

    corr.main()
