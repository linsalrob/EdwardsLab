import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import cartopy.crs as ccrs
import re

from roblib import parse_dnadist
from roblib import latlon2distance

import sys

"""
Read the dna distances (note file names are fixed at the moment), latitude, and longitude of the samples
and plot a line based on that
"""

def get_lon_lat(idf, maxtoget=50000):
    """
    Get the longitude and latitude of different ids. Note that we have longitude first to work with cartopy
    :param idf: the id.map file
    :param maxtoget: the maxiumum number of ids to get. This is just for debugging
    :return:
    """
    lonlat = {}
    count = 0
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
            lonlat[p[0]] = (lon, lat)
    return lonlat

def get_dna_distance(ddf):
    """
    Get the DNA distance between pairs of IDs
    :param ddf: dna distance file (e.g. from neighbor)
    :return: hash of {id1}->{id2}->{distance}
    """

    ids, distance = parse_dnadist(ddf)

    ids = [x.replace('_R_', '') for x in ids]

    # we want similarity (ie. close = 1) not distance
    # to highlight closer things
    sims = []
    for i in range(len(distance)):
        sims.append([1-float(x) for x in distance[i]])

    similarities = {}
    for idx1, name1 in enumerate(ids):
        similarities[name1] = {}
        for idx2, name2 in enumerate(ids):
            similarities[name1][name2] = sims[idx1][idx2]
    return similarities


def best_dna_dist(ddf, kms={}):
    """
    Get the best dna distance between two samples.

    :param ddf: dnadist file
    :param kms: a hash of hashes of distance (e.g. in kms) between the samples with IDs. If possible we use this to
                resolve duplicately lowest scores
    :return:
    """

    ids, distance = parse_dnadist(ddf)
    ids = [x.replace('_R_', '') for x in ids]

    bestdists = {}
    for idx, name in enumerate(ids):
        # of course  the minimum distance is itself and we need to ignore that
        mindist = 100
        bestdists[name]={}
        bestname = None
        for posn, dist in enumerate([float(x) for x in distance[idx]]):
            if name not in kms or ids[posn] not in kms[name]:
                continue
            if kms[name][ids[posn]] == 0:
                # we are at the same site
                continue
            if dist < mindist and posn != idx:
                bestname = ids[posn]
                mindist = dist
            if dist == mindist and posn != idx and bestname and \
                            name in kms and bestname in kms[name] \
                            and kms[name][ids[posn]] < kms[name][bestname]:
                bestname = ids[posn]

        bestdists[name][bestname] = 1 - mindist

    return bestdists



def plotmap(ll, dd, maxdist=1, maxlinewidth=3):
    """
    Plot the map of the dna distances and lat longs
    :param ll: The lon-lats
    :param dd: The distances to use
    :param maxdist: The maximum distance that we will scale to be maxlinewidth
    :return:
    """

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
    samples = ll.keys()

    # plot the circles for each sample site
    for lonlat in ll.values():
        plt.plot(lonlat[0], lonlat[1], '*', markersize=18, transform=ccrs.PlateCarree())

    for idx1 in samples:
        for idx2 in samples:
            if idx1 == idx2:
                continue
            # this should only happen when we do best DNA distances
            if idx2 not in dd[idx1]:
                continue
            linewidth = dd[idx1][idx2]
            linewidth = linewidth/maxdist * maxlinewidth
            colorVal = scalarMap.to_rgba(dd[idx1][idx2])
            plt.plot([ll[idx1][0], ll[idx2][0]], [ll[idx1][1], ll[idx2][1]], color=colorVal, linewidth=linewidth, alpha=0.25, transform=ccrs.Geodetic())

#    plt.colorbar(CS3)

    plt.show()


if __name__ == '__main__':
    dnadist_file = "/home/redwards/Dropbox/GitHubs/crAssphage/Global_Survey/Analysis/PrimerA/seqs.dnadist"
    idmap_file = "/home/redwards/Dropbox/GitHubs/crAssphage/Global_Survey/Analysis/PrimerA/id.map"
    lonlat = get_lon_lat(idmap_file)
    # This line gets _all_ similarites. The next line gets the single best
    # similarities = get_dna_distance(dnadist_file)

    # convert the latlon to km
    kms = {}
    for name1 in lonlat.keys():
        kms[name1] = {}
        for name2 in lonlat.keys():
            # be careful here about the order of lat and lon
            kms[name1][name2] = latlon2distance(float(lonlat[name1][1]), float(lonlat[name1][0]),
                                                float(lonlat[name2][1]), float(lonlat[name2][0]))
    similarities = best_dna_dist(dnadist_file, kms)
    plotmap(lonlat, similarities)
