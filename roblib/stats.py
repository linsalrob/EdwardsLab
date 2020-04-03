import numpy
import sys
from math import sqrt

__author__ = 'Rob Edwards'



def mean(arr):
    """
    Calculate the mean of all the values in a list

    :param arr: The list of values
    :type arr: list
    :return: The mean of all the numbers
    :rtype: float
    """
    n = len(arr)
    sm = sum(arr)
    return float(sm) / float(n)


def median(arr):
    """
    Calculate the median of all the values in a list

    :param arr: The list of values
    :type arr: list
    :return: The median of all the numbers
    :rtype: float
    """
    arr.sort()
    n = len(arr)
    mid = int(n / 2)
    return arr[mid]


def stdev(arr):
    """
    Calculate the standard deviation of all the values in a list

    :param arr: The list of values
    :type arr: list
    :return: The st dev of all the numbers
    :rtype: float
    """
    n = len(arr)
    sm = sum(arr)
    nmpyarr = numpy.array(arr)
    sum_of_squares = numpy.sum(nmpyarr * nmpyarr)
    result = 0
    try:
        result = sqrt(sum_of_squares / n - (sm / n) ** 2)
    except:
        sys.stderr.write("can't create the stdev. Only have " + str(n) + " elements in the array")

    return result

