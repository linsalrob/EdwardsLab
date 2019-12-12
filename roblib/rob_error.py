"""
An Error Class so I can write my own errors
"""
class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass

class SequencePairError(Error):
    """
    Exception raised for sequences not being paired properly.

    :param message: explanation of the error
    """

    def __init__(self, message):
        self.message = message

class FastqFormatError(Error):
    """
    Exception raised for sequences not being paired properly.

    :param message: explanation of the error
    """

    def __init__(self, message):
        self.message = message
