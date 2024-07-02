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

    def __init__(self, input_cn, input_gc_map, output, mappability=0.9,
                 smoothing_function='lowess',
                 polynomial_degree=2):
        self.mappability = mappability

        self.input_cn = input_cn
        self.input_gc_map = input_gc_map
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
        # read input dataframe containing simulated read counts and copy number for each bin
        input_cn = pd.read_csv(self.input_cn, sep='\t')

        # rename true_reads_raw to reads and subset to the input columns for hmmcopy
        input_cn.rename(columns={'true_reads_raw': 'reads'}, inplace=True)
        input_cn = input_cn[['chr', 'start', 'end', 'reads', 'cell_id']]
        print('input_cn', input_cn.head(), input_cn.dtypes, sep='\n')

        # read input for gc and mappability of each bin
        input_gc_map = pd.read_csv(self.input_gc_map)
        # subtract 1 from the start position so the bins match input_cn
        input_gc_map['start'] = input_gc_map['start'] - 1
        print('input_gc_map', input_gc_map.head(), input_gc_map.dtypes, sep='\n')

        # merge input_cn and input_gc_map to create input_df
        # this dataframe has columns 'chr', 'start', 'end', 'gc', 'map', 'reads', 'cell_id'
        input_df = pd.merge(input_cn, input_gc_map)
        print(input_df.head())

        # store each cell's dataframe in a list of dataframes
        output_df = []

        # iterate over each cell
        for cell_id, cell_df in input_df.groupby('cell_id'):
            print("Processing cell: {}".format(cell_id))
            print(cell_df.head())

            # create a copy of the cell dataframe with reset index
            df = cell_df.copy().reset_index(drop=True)

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

            print('done with cell: {}'.format(cell_id))
            print(df.head())

            output_df.append(df)
    
        # concatenate the per-cell dataframes into one output dataframe
        print('number of cells processed: {}'.format(len(output_df)))
        output_df = pd.concat(output_df, ignore_index=True)

        # save
        self.write(output_df)


def parse_args():
    """
    parses command line arguments
    """

    parser = argparse.ArgumentParser()

    parser.add_argument('input_cn',
                        help='path to the input read count csv'
                        )

    parser.add_argument('input_gc_map',
                        help='path to the input gc and mappability csv'
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

    corr = CorrectReadCount(args.input_cn, args.input_gc_map, args.output,
                            mappability=args.mappability,
                            )

    corr.main()
