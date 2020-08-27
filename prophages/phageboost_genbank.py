"""
Run phage boost.

This is extracted from their example code, but we added parsing genbank files to dataframes

"""

import os
import sys
from roblib import genbank_to_pandas, message
import pandas as pd
from PhageBoost.main import read_sequence_file, call_genes, calculate_features, read_model_from_file, predict, get_predictions
import xgboost as xgb
import argparse

__author__ = 'Rob Edwards'
__copyright__ = 'Copyright 2020, Rob Edwards'
__credits__ = ['Rob Edwards']
__license__ = 'MIT'
__maintainer__ = 'Rob Edwards'
__email__ = 'raedwards@gmail.com'


def run_phage_boost(genecalls, model_file, verbose):
    """
    Run phage boost
    :param model_file: The model file that is probably something like model_delta_std_hacked.pickled.silent.gz
    :param genecalls: The pandas data frame of gene calls
    :param verbose: more output
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

    # calculate features from gene calls
    if verbose:
        message("Calculating features", "GREEN")

    df = calculate_features(genecalls)
    # load model
    model, feats, feats_, limit = read_model_from_file(model_file)
    # transform data
    df = get_predictions.get_deltas(df[feats_])
    if verbose:
        message("Transforming gene predictions to regions", "GREEN")
    # transform single gene predictions to regions
    newgenecalls, nphages, res = predict(model, genecalls, df,
                                         feats, period, win_type,
                                         min_periods, limit, threshold,
                                         length, gaps, neighbouring, alpha)
    return res



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=" ")
    parser.add_argument('-g', '--genbankfile', help='GenBank file to parse', required=True)
    parser.add_argument('-m', '--modelfile', required=True,
                        help="Model file. Probably something like model_delta_std_hacked.pickled.silent.gz")
    parser.add_argument('-o', '--outputfile', help='output file for phage regions')
    parser.add_argument('-c', '--mincontiglen', default=1000, type=int,
                        help='minimum contig length  [Default: %(default)d]')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    if args.verbose:
        message("Reading genbank file", "GREEN")
    genecalls = genbank_to_pandas(args.genbankfile, args.mincontiglen, True, True, args.verbose)
    if args.verbose:
        message("Phage Boosting", "GREEN")
    res = run_phage_boost(genecalls, args.modelfile, args.verbose)
    if args.outputfile:
        with open(args.outputfile, 'w') as out:
            res.to_csv(out, sep="\t", header=True)
    else:
        print(res)
