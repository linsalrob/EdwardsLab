"""
Run phage boost.

This is extracted from their example code, but we added parsing genbank files to dataframes

"""

import os
import sys
from roblib import genbank_to_pandas
from PhageBoost.main import read_sequence_file, call_genes, calculate_features, read_model_from_file, predict, get_predictions
import xgboost as xgb
import argparse

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'

def run_phage_boost(genecalls, verbose):
    """
    Run phage boost
    :param genecalls:
    :param verbose:
    :return:
    """
    # rolling params
    period = 20
    win_type = 'parzen'
    min_periods = 1

    # region finding params
    threshold = 0.9
    length = 10
    gaps = 5
    neighbouring = 0
    alpha = 0.001

    # files
    model_file = '/home/redwards/pb/PhageBoost/PhageBoost/models/model_delta_std_hacked.pickled.silent.gz'

    # calculate features from gene calls
    df = calculate_features(genecalls)
    # load model
    model, feats, feats_, limit = read_model_from_file(model_file)
    # transform data
    df = get_predictions.get_deltas(df[feats_])
    # transform single gene predictions to regions
    newgenecalls, nphages, res = predict(model, genecalls, df,
                                         feats, period, win_type,
                                         min_periods, limit, threshold,
                                         length, gaps, neighbouring, alpha)
    return res



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', help='GenBank file to parse', required=True)
    parser.add_argument('-o', help='output file for phage regions')
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    genecalls = genbank_to_pandas(args.g, args.v)
    res = run_phage_boost(genecalls, args.v)
    if args.o:
        with open(args.o, 'w') as out:
            res.to_csv(out, sep="\t", header=True)
    else:
        print(res)