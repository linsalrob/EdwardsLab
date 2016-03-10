"""
Resample 80% of the data and plot a graph of how many new things we see. This is to answer an argument with Geni

"""

import os
import sys
import argparse
import matplotlib.pyplot as plt
from random import shuffle

def resample(size, percent, tries):
    if percent > 1:
        percent /= 100

    # define an array of size size
    data = [i for i in range(size)]
    # where we put the results as a cumulative total
    iterations = []
    seen = set()

    for t in range(tries):
        # randomize the array
        shuffle(data)
        # see if we have seen percent things
        new = 0
        resampsize = int(size * percent)
        # sys.stderr.write("resampling " + str(resampsize) + " from " + str(size) + "\n")
        for i in range(resampsize):
            if data[i] not in seen:
                seen.add(data[i])
                new += 1

        if not iterations:
            iterations.append(new)
        else:
            iterations.append(new+iterations[-1])

    # now just plot the number of new things as a cumulative total
    plt.plot(iterations)
    plt.ylabel('New numbers seen')
    plt.xlabel('Iteration')
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Resample a list of numbers to see the new things seen")
    parser.add_argument('-s', help='Size of array to resample from (size of dataset)', type=int, required=True)
    parser.add_argument('-p', help='Percent to resample at each iteration (float)', type=float, required=True)
    parser.add_argument('-i', help='Number of iterations to run', type=int, required=True)
    args = parser.parse_args()

    resample(args.s, args.p, args.i)