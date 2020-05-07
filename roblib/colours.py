"""
Colors that you can import and make the text look pretty

Source: https://stackoverflow.com/questions/287871/print-in-terminal-with-colors
"""

import sys

__author__ = 'Rob Edwards'


class colours(object):
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    PINK = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE = '\033[0m'


class colors(object):
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

    PINK = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    WHITE = '\033[0m'


def message(msg, color):
    """
    Print a message to stderr using color
    :param msg: the message to print
    :param color: the color to use
    :return: nothing
    """

    colours = {
        "BLUE" : '\033[94m',
        "BOLD" : '\033[1m',
        "ENDC" : '\033[0m',
        "FAIL" : '\033[91m',
        "GREEN" : '\033[92m',
        "HEADER" : '\033[95m',
        "PINK" : '\033[95m',
        "RED" : '\033[91m',
        "UNDERLINE" : '\033[4m',
        "WARNING" : '\033[93m',
        "WHITE" : '\033[0m',
        "YELLOW" : '\033[93m',
    }

    color = color.upper()
    if color not in colours:
        colours[color] = '\033[0m'
        sys.stderr.write(f"Error: No colour {color}\n")

    sys.stderr.write(f"{colours[color]}{msg}{colours['ENDC']}\n")