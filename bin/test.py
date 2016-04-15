"""
Test stuff!
"""

import os
import sys

import argparse
import scipy
import scipy.cluster.hierarchy as sch

import dateutil.parser

d = 'Thu Feb 11 16:39:56 -0800 2016'
z = dateutil.parser.parse(d)
print(z)




sys.exit(0)


X = scipy.randn(100, 2)  # 100 2-dimensional observations
print(X)

d = sch.distance.pdist(X)
print(len(d))
