import matplotlib.pyplot as plt

import cartopy.crs as ccrs

"""
This example is from the cartopy website, and is mostly to make sure things are installed and working.
"""


def main():
    ax = plt.axes(projection=ccrs.Robinson())

    # make the map global rather than have it zoom in to
    # the extents of any plotted data
    ax.set_global()

    ax.stock_img()
    ax.coastlines()

    # san diego
    sdlat, sdlon = 32.7157, -117.1611
    # brisbane
    brislat, brislon = -27.4698, 153.0251


    # NOTE: longitude before latitude!!
    plt.plot([sdlon, brislon], [sdlat, brislat], color='blue', linewidth=2,  transform=ccrs.Geodetic())



    plt.show()


if __name__ == '__main__':
    main()
