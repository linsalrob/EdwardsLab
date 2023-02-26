"""
Some settings for the config files
"""

import os
import sys


def get_db_dir():
    """
    Just return the default dir listed above
    :return: the default location for the sqllite database
    """

    if 'NCBI_TAXONOMY' in os.environ:
        if os.path.exists(os.environ['NCBI_TAXONOMY']):
            return os.environ['NCBI_TAXONOMY']
        else:
            print(f"WARNING: NCBI_TAXONOMY variable is set but {os.environ['NCBI_TAXONOMY']} does not exist", file=sys.stderr)
    if 'TAXONKIT_DB' in os.environ:
        if os.path.exists(os.environ['TAXONKIT_DB']):
            return os.environ['TAXONKIT_DB']
        else:
            print(f"WARNING: TAXONKIT_DB variable is set but {os.environ['TAXONKIT_DB']} does not exist", file=sys.stderr)

    return None
