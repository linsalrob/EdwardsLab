"""
Some settings for the config files
"""

# defaultdir = '/data/ncbi/taxonomy/current'
# defaultdir = '/home/edwa0468/ncbi/taxonomy'
defaultdir = '/raid60/usr/data/NCBI/taxonomy/current/'

def get_db_dir():
    """
    Just return the default dir listed above
    :return: the default location for the sqllite database
    """
    return defaultdir
