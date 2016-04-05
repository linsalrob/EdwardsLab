"""
Test stuff!
"""

import os
import sys

import argparse
import scipy
import scipy.cluster.hierarchy as sch


X = scipy.randn(100, 2)  # 100 2-dimensional observations
print(X)

d = sch.distance.pdist(X)
print(len(d))
