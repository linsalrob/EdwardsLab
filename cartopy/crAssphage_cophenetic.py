"""
Read a cophenetic matrix to plot the distances. You can make the matrix using ete3 and tree_to_cophenetic_matrix.py
"""

import os
import sys
import argparse

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import math

import cartopy.crs as ccrs
import re

def get_lon_lat(idf, maxtoget=50000):
    """
    Get the longitude and latitude of different ids. Note that we have longitude first to work with cartopy
    :param idf: the id.map file
    :param maxtoget: the maxiumum number of ids to get. This is just for debugging
    :return:
    """
    lonlat = {}
    count = 0
    global verbose
    with open(idf, 'r') as fin:
        for l in fin:
            if count > maxtoget:
                break
            count+=1
            s=re.search('latitude=(\S+)\]', l)
            if not s:
                sys.stderr.write("No latitude in {}".format(l))
                continue
            lat=s.group(1)

            s = re.search('longitude=(\S+)\]', l)
            if not s:
                sys.stderr.write("No longitude in {}".format(l))
                continue
            lon = s.group(1)
            p=l.split("\t")
            lonlat[p[0]] = (float(lon), float(lat))
    return lonlat

def latlon2distance(lat1, long1, lat2, long2, miles=False):
    """Convert two coordinates to distance.

    This is an approximation since the earth is not spherical, but accuracy is <100m, especially for close points

    This code was taken from http://www.johndcook.com/python_longitude_latitude.html

    Latitude is measured in degrees north of the equator; southern locations have negative latitude.
    Similarly, longitude is measured in degrees east of the Prime Meridian. A location 10deg west of
    the Prime Meridian, for example, could be expressed as either 350deg  east or as -10deg east.

    Arguments: lat1, long1; lat2, long2; miles is a boolean. If you want miles set it to true. Else set it to false

    """
    global verbose

    if lat1 == lat2 and long1 == long2:
        return 0


    # Convert latitude and longitude to
    # spherical coordinates in radians.
    degrees_to_radians = math.pi / 180.0

    # phi = 90 - latitude
    phi1 = (90.0 - lat1) * degrees_to_radians
    phi2 = (90.0 - lat2) * degrees_to_radians

    # theta = longitude
    theta1 = long1 * degrees_to_radians
    theta2 = long2 * degrees_to_radians

    # Compute spherical distance from spherical coordinates.

    # For two locations in spherical coordinates
    # (1, theta, phi) and (1, theta, phi)
    # cosine( arc length ) =
    #    sin phi sin phi' cos(theta-theta') + cos phi cos phi'
    # distance = rho * arc length

    cos = (math.sin(phi1) * math.sin(phi2) * math.cos(theta1 - theta2) + math.cos(phi1) * math.cos(phi2))
    try:
        arc = math.acos(cos)
    except Exception as err:
        sys.stderr.write("There was an err: {} trying to take the acos of ({})\n".format(err, cos))
        arc=0
    # Remember to multiply arc by the radius of the earth
    # in your favorite set of units to get length.
    #
    # To convert to miles multiple arc by 3960
    # To convert to kilometers multiply arc by 6373

    if miles:
        arc *= 3960
    else:
        arc *= 6373

    return arc


def closest_dna_dist(matrixfile):
    """
    Read the matrix file and get the id of the point with the closest distance that is not ourself
    :param treefile: The cophenetic matrix file to read
    :return: a dict of a node and its closest leaf
    """

    global verbose
    if verbose:
        sys.stderr.write("Getting closest distances\n")
    distances = {}
    with open(matrixfile, 'r') as f:
        l = f.readline()
        ids = l.rstrip().split("\t")
        for i,name in enumerate(ids):
            if i == 0:
                continue
            distances[name] = {}
        for l in f:
            data = l.rstrip().split("\t")
            for i,dist in enumerate(data):
                if i == 0:
                    continue
                distances[data[0]][ids[i]] = float(dist)
                distances[ids[i]][data[0]] = float(dist)

    closest = {}
    for d in distances:
        closest[d] = {}
        for k in sorted(distances[d], key=distances[d].get):
            if k == d:
                continue
            closest[d][k] = distances[d][k]
            break

    if verbose:
        sys.stderr.write("Done\n")
    return closest

def plotmap(ll, dd, outputfile, maxdist=1, maxlinewidth=3):
    """
    Plot the map of the dna distances and lat longs
    :param ll: The lon-lats
    :param dd: The distances to use
    :param outputfile: The file name to write the image to
    :param maxdist: The maximum distance that we will scale to be maxlinewidth
    :return:
    """
    global verbose

    ax = plt.axes(projection=ccrs.Robinson())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()

    ax.stock_img()
    ax.coastlines()

    ## color the lines based on the maximum distance value
    jet = cm = plt.get_cmap('jet')
    cNorm = colors.Normalize(vmin=0, vmax=maxdist)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

    # Using contourf to provide my colorbar info, then clearing the figure
    Z = [[0, 0], [0, 0]]
    levels = range(0, int(100 * maxdist) + 10, 10)
    CS3 = plt.contourf(Z, levels, cmap=jet)
#    plt.clf()




    # NOTE: longitude before latitude!!
    # plt.plot([sdlon, brislon], [sdlat, brislat], color='blue', linewidth=2,  transform=ccrs.Geodetic())

    # plot the circles for each sample site
    # markerfacecolor="None",
    for lonlat in ll.values():
        plt.plot(lonlat[0], lonlat[1], 'o', color='Black', alpha=0.25, markersize=4, transform=ccrs.PlateCarree())

    for idx1 in dd:
        for idx2 in dd[idx1]:
            # this should only happen when we do best DNA distances
            if idx1 not in ll:
                sys.stderr.write("NO Lat/Lon for {}\n".format(idx1))
                continue
            if idx2 not in ll:
                sys.stderr.write("NO Lat/Lon for {}\n".format(idx2))
                continue

            if verbose:
                sys.stderr.write("Distance between {} and {}: {}\n".format(idx1, idx2, latlon2distance(ll[idx1][1], ll[idx1][0], ll[idx2][1], ll[idx2][0])))

            linewidth = dd[idx1][idx2]
            linewidth = linewidth/maxdist * maxlinewidth
            #colorVal = scalarMap.to_rgba(dd[idx1][idx2])
            plt.plot([ll[idx1][0], ll[idx2][0]], [ll[idx1][1], ll[idx2][1]], color='Red', linewidth=linewidth, alpha=0.1, transform=ccrs.Geodetic())

            if latlon2distance(ll[idx1][1], ll[idx1][0], ll[idx2][1], ll[idx2][0]) < 100:
                if verbose:
                    sys.stderr.write("Adding a circle for {} and {}\n".format(ll[idx1][0], ll[idx1][1]))
                plt.plot(ll[idx1][0], ll[idx1][1], 'o', color='Red', alpha=0.1, markersize=2,
                         transform=ccrs.PlateCarree())




#    plt.colorbar(CS3)

    #plt.show()
    plt.savefig(outputfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot a map using ete and lat lon')
    parser.add_argument('-i', help='id.map file with lat/lon information', required=True)
    parser.add_argument('-m', help='cophenetic map file with same ids as id.map', required=True)
    parser.add_argument('-o', help='output file name', required=True)
    parser.add_argument('-v', help='verbose output', action='store_true')
    args = parser.parse_args()

    global verbose
    verbose = False
    if args.v:
        verbose = True

    lonlat = get_lon_lat(args.i)
    # dist = best_dna_dist(get_dna_distance(args.t))
    dist = closest_dna_dist(args.m)
    plotmap(lonlat, dist, args.o)