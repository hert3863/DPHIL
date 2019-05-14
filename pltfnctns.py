# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 13:46:51 2016

@author: bakerh

NAME
    Plotting functions
PURPOSE
    Defines all plotting functions for AMO
"""


import numpy as np


def timeplot(series, startdate, mnths, colour):
    '''
    timeplot takes a times series and outputs a plot of the timeseries from
    the start date

    Parameters
    ----------
    series : array
        Time series to be plotted
    startdate : int
        time series start date
    mnths : int
        number of months per year in data
    colour : str
        colour of data line

    Returns
    -------
    plot : fig
        plot of time series
    '''
    import matplotlib.pyplot as plt
    import datetime as dt
    plt.figure()
    ts = np.array([dt.date(i, j, 1) for i in range(startdate, startdate +
                  int(len(series)/mnths)) for j in range(1, mnths + 1)])
    plt.plot(ts, series, color=colour)
    plt.axhline(0, color='black')
    plt.show()


def lanczosfilter(timeseries, window, cutoff):
    """
    Calculate weights for a low pass Lanczos filter
    and then applies to time series

    Parameters
    ----------
    timeseries: array
        Time series to be filtered
    window: int
        The length of the filter window.
    cutoff: float
        The cutoff frequency in inverse time steps.

    Returns
    -------
    filteredtimeseries: array
        low pass filtered time series
    """
    order = ((window - 1) // 2) + 1
    nwts = 2 * order + 1
    w = np.zeros([nwts])
    n = nwts // 2
    w[n] = 2 * cutoff
    k = np.arange(1., n)
    sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
    firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
    w[n-1:0:-1] = firstfactor * sigma
    w[n+1:-1] = firstfactor * sigma
    w = w[1:-1]
    # Set up loop to run across all time data
    filteredtimeseries = np.array([0])
    for i in range(len(timeseries)-len(w)+1):
        # Set up loop to compute filter*data
        weightsum = 0
        for j in range(i, i+len(w)):
            weighted = timeseries[j] * w[j-i]
            weightsum = weightsum + weighted
        ws = np.array([weightsum])
        filteredtimeseries = np.concatenate((filteredtimeseries, ws))
    return filteredtimeseries[1:]


def corrmap(timeseries, ts, griddeddata, gs, mnths):
    """
    Produces a correlation grid between time series and gridded data
    Takes a time series and correlates it with every grid point time
    series, outputting a correlation grid and significance grid. Time
    series and gridded data should have the same length of time

    Parameters
    ----------
    timeseries: array
        Time series to be correlated
    ts: int
        year time series starts
    griddeddata: array
        Gridded data to be correlated
    gs: int
        year gridded data starts
    mnths: int
        number of months per year in both time series

    Returns
    -------
    corrmap: array
        map of correlated data
    sigmap: array
        mag of significance of correlated data
    """
    import scipy.stats
    # make two datasets same time length
    if ts < gs:
        timeseries = timeseries[(gs-ts)*mnths:-1]
    elif ts > gs:
        griddeddata = griddeddata[(ts-gs)*mnths:-1]
    if np.ma.size(griddeddata, axis=0) < len(timeseries):
        timeseries = timeseries[0:(np.ma.size(griddeddata, axis=0))]
    elif np.ma.size(griddeddata, axis=0) > len(timeseries):
        griddeddata = griddeddata[0:len(timeseries), :, :]
    # calculate required grids
    corr = np.zeros([np.ma.size(griddeddata, axis=1),
                     np.ma.size(griddeddata, axis=2)])
    sig = np.zeros([np.ma.size(griddeddata, axis=1),
                    np.ma.size(griddeddata, axis=2)])
    for i in range(np.ma.size(griddeddata, axis=1)):
        for j in range(np.ma.size(griddeddata, axis=2)):
            (c, s) = scipy.stats.pearsonr(timeseries, griddeddata[:, i, j])
            corr[i, j] = c
            sig[i, j] = s
    return (corr, sig)


def seasonseries(timeseries):
    """
    Takes a time series and returns two timeseries, one
    of the winter months (DJF) and one of the summer
    months (JJA)

    Parameters
    ----------
    timeseries: array
        timeseries to be split

    Returns
    -------
    winter: array
        DJF time series
    summer: array
        JJA time series
    """
    winter = np.zeros([int(len(timeseries)/4)-3])
    summer = np.zeros([int(len(timeseries)/4)])
    for i in range(1, len(timeseries)-11, 12):
        for j in range(11, 14):
            winter[int((i-1)/4)+j-11] = timeseries[i+j-1]
    for i in range(1, len(timeseries)+1, 12):
        for j in range(5, 8):
            summer[int((i-1)/4)+j-5] = timeseries[i+j-1]
    return winter, summer


def gridseasonseries(grid):
    """
    Takes a gridded time series and returns two griddedtimeseries, one
    of the winter months (DJF) and one of the summer
    months (JJA)

    Parameters
    ----------
    grid: array
        gridded data time series to be splot into seasons

    Returns
    -------
    winter: array
        DJF time series
    summer: array
        JJA time series
    """
    winter = np.zeros([int(np.ma.size(grid, axis=0)/4)-3,
                       np.ma.size(grid, axis=1), np.ma.size(grid, axis=2)])
    summer = np.zeros([int(np.ma.size(grid, axis=0)/4),
                       np.ma.size(grid, axis=1), np.ma.size(grid, axis=2)])
    for a in range(np.ma.size(grid, axis=1)):
        for b in range(np.ma.size(grid, axis=2)):
            for i in range(1, (np.ma.size(grid, axis=0)-11), 12):
                for j in range(11, 14):
                    winter[int((i-1)/4)+j-11, a, b] = grid[i+j-1, a, b]
    for a in range(np.ma.size(grid, axis=1)):
        for b in range(np.ma.size(grid, axis=2)):
            for i in range(1, np.ma.size(grid, axis=0)+1, 12):
                for j in range(5, 8):
                    summer[int((i-1)/4)+j-5, a, b] = grid[i+j-1, a, b]
    return winter, summer


def colourscale(plotdata):
    """
    Takes data being plotted and normalises the colourscale between largest
    data value and it's negative multiple

    Parameters
    ----------
    plotdata: array
        data being plotted

    Returns
    -------
    caxismax: int
        max magnitude data value
    caxismin: int
        negative of max mag data value
    ctrs: array
        gradated colour scale contour array
    """
    M = np.nanmax(plotdata)
    m = np.nanmin(plotdata)
    if M >= abs(m):
        ctrs1 = np.arange(-M, 0, .1*M)
        ctrs2 = np.arange(0.1*M, 1.1*M, .1*M)
        ctrs = np.concatenate((ctrs1, ctrs2))
        caxismin = -M
        caxismax = M
    else:
        m = -m
        ctrs1 = np.arange(-m, 0, .1*m)
        ctrs2 = np.arange(0.1*m, 1.1*m, .1*m)
        ctrs = np.concatenate((ctrs1, ctrs2))
        caxismin = -m
        caxismax = m
    # function will not work if there exist no positive max or negative min
    return caxismin, caxismax, ctrs


def mapplot(plotdata, title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap, shiftgrid
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    lat = np.arange(90, -90-180/(np.ma.size(plotdata, axis=0)-1),
                    -180/(np.ma.size(plotdata, axis=0)-1))
    lon = np.arange(0, 360+360/(np.ma.size(plotdata, axis=1)),
                    360/(np.ma.size(plotdata, axis=1)))
    temp = plotdata[:, 0]
    plotdata = np.c_[plotdata, temp]
    plotdata, lon = shiftgrid(180., plotdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)
    plt.figure()
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    my_cmap = plt.cm.jet(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscale(plotdata)
    plot = m.contourf(x, y, plotdata, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def mapplotsst(plotdata, title=''):
    """
    Plots input Kaplan SST grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        SST data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap, shiftgrid
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    lat = np.arange(-87.5, 87.5+180/(np.ma.size(plotdata, axis=0)),
                    180/(np.ma.size(plotdata, axis=0)))
    lon = np.arange(2.5, 362.5+360/(np.ma.size(plotdata, axis=1)),
                    360/(np.ma.size(plotdata, axis=1)))
    temp = plotdata[:, 0]
    plotdata = np.c_[plotdata, temp]
    plotdata, lon = shiftgrid(180., plotdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)
    plt.figure()
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    my_cmap = plt.cm.jet(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscale(plotdata)
    plotdata = np.ma.masked_invalid(plotdata, copy=True)
    plt.gca().patch.set_color('.75')
    plot = m.contourf(x, y, plotdata, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def jetindices(timeseriesgrid):
    """
    Outputs maximum uwnd speed and position of maximum over
    North Atlantic region (20-70N and 0-90W) for a gridded time series
    Inputted gridded data over many months.

    MUST BE 20thC RENALYSIS FOR DIMENSIONS TO WORK

    Parameters
    ----------
    timeseriesgrid: array
        data to determine jet indices for

    Returns
    -------
    jetspeed: array
        max jet speed
    jetlatitude: array
        position of maximum
    """
    temp = timeseriesgrid[:, :, 0:1]
    timeseriesgrid = np.c_[timeseriesgrid, temp]
    atlanticgrid = timeseriesgrid[:, 10:36, 135:]

    def indices(grid):
        m = np.amax(grid, 1)
        speed = np.amax(m)
        ind = np.argmax(m, 0)
        lat = 70-ind*2
        return speed, lat
    positions = np.zeros(1)
    speeds = np.zeros(1)
    for i in range(np.ma.size(atlanticgrid, axis=0)):
        speed, lat = indices(atlanticgrid[i, :, :])
        s, l = np.array([speed]), np.array([lat])
        positions = np.concatenate((positions, l))
        speeds = np.concatenate((speeds, s))
    positions = positions[1:]
    speeds = speeds[1:]
    return positions, speeds
