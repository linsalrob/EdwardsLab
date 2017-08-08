"""

"""

import os
import sys
import argparse
import bz2
from bs4 import BeautifulSoup



# with bz2.BZ2File('biosample_set.xml.bz2', 'r') as input:
#     for i in range(10):
#         l = input.readline()
#         print("{}\n".format(l))


# <Id db="BioSample" is_primary="1">SAMN00000002</Id>\n'

def primaryId(tag):
    return tag['db'] == 'BioSample' and tag['is_primary']


with bz2.BZ2File('biosample_set.xml.bz2', 'r') as input:
    soup = BeautifulSoup(input, 'xml')
    pi = soup.find_next(primaryId)
    print("{}".format(pi))
