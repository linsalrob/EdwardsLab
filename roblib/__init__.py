import os
import sys

__author__ = 'Rob Edwards'
from .stats import mean, median, stdev
from .sequences import read_fasta, readFasta, stream_fastq
from .dna import rc, shannon
from .geography import latlon2distance
from .strings import ascii_clean
from .functions import is_hypothetical
from .newick import Newick_Tree

__all__ = [
    'mean', 'median', 'stdev',
    'read_fasta', 'readFasta', 'stream_fastq',
    'rc', 'shannon',
    'latlon2distance',
    'ascii_clean', 'is_hypothetical', 'Newick_Tree'
    ]
