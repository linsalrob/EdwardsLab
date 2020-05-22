"""
Custom exceptions for taxnomy parsing
"""

class Error(Exception):
    """Base class for other exceptions"""
    pass


class EntryNotInDatabaseError(Exception):
    """Entry not in the db. Obvs"""

    def __init__(self, message):
        self.message = message


class NoNameFoundError(Exception):
    """No name was found for this entry"""
    def __init__(self, message):
        self.message = message