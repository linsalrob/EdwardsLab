import os
import sys

__author__ = 'Rob Edwards'
from .stats import mean, median, stdev
from .sequences import read_fasta, readFasta, stream_fastq, stream_fasta, stream_paired_fastq
from .dna import rc, shannon
from .geography import latlon2distance
from .strings import ascii_clean
from .functions import is_hypothetical
from .newick import Newick_Tree
from .dnadist import parse_dnadist
from .blast import stream_blast_results
from .translate import translate_dna
from .bcolors import bcolors
from .rob_error import SequencePairError, FastqFormatError

__all__ = [
    'mean', 'median', 'stdev',
    'read_fasta', 'readFasta', 'stream_fastq', 'stream_fasta', 'stream_paired_fastq',
    'rc', 'shannon',
    'latlon2distance',
    'ascii_clean', 'is_hypothetical', 'Newick_Tree', 'parse_dnadist',
    'stream_blast_results',
    'translate_dna',
    'bcolors',
    'SequencePairError', 'FastqFormatError'
    ]
