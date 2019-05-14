# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 16:18:01 2017

Klingaman analysis

@author: bakerh
"""

import numpy as np


def ncread(filelocation, invariable):
    '''
    ncread outputs numpy arrays of called variables in netCDF4 file

    Parameters
    ----------
    fileloaction : str
        NetCDFfile
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variable : array
        Array of variable specified
    '''
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    nc_f = filelocation  # filename
    nc_fid = Dataset(nc_f, 'r')
    variable = nc_fid.variables[invariable][:]
    return variable


def reorder(inp):
    out = np.zeros((np.ma.size(inp, axis=0),
                    np.ma.size(inp, axis=1), np.ma.size(inp, axis=2)))
    for i in range(np.floor_divide(np.ma.size(inp, axis=0), 12)):
        out[1+i*12, :] = inp[4+i*12, :]
        out[2+i*12, :] = inp[3+i*12, :]
        out[3+i*12, :] = inp[7+i*12, :]
        out[4+i*12, :] = inp[0+i*12, :]
        out[5+i*12, :] = inp[8+i*12, :]
        out[6+i*12, :] = inp[6+i*12, :]
        out[7+i*12, :] = inp[5+i*12, :]
        out[8+i*12, :] = inp[1+i*12, :]
        out[9+i*12, :] = inp[11+i*12, :]
        out[10+i*12, :] = inp[10+i*12, :]
        out[11+i*12, :] = inp[9+i*12, :]
        out[0+i*12, :] = inp[2+i*12, :]
    return out


def processall():
    import os
    from netCDF4 import Dataset
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')
    lsm = lsm[0, 0, :]
    control = np.zeros((4, 144, 192))
    for d in os.listdir('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/'):
        if d.find('control') == 0:
            nc_f = d
            nc_fid = Dataset('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/' + nc_f, 'r')
            cont = nc_fid.variables['air_temperature'][480:3000, :]
            cont = np.squeeze(cont)
            cont = reorder(cont)
            control[0, :] += np.mean(cont[0:3, :], axis=0)
            control[1, :] += np.mean(cont[3:6, :], axis=0)
            control[2, :] += np.mean(cont[6:9, :], axis=0)
            control[3, :] += np.mean(cont[9:12, :], axis=0)
            print('Done ' + d)
    control = control/10
    for k in range(4):
        control[k, :] = np.ma.masked_array(control[k, :],
                                           np.logical_not(lsm[144:0:-1, :]))

    zero98 = np.zeros((4, 144, 192))
    for d in os.listdir('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/'):
        if d.find('0p98.') == 0:
            nc_f = d
            nc_fid = Dataset('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/' + nc_f, 'r')
            zero = nc_fid.variables['air_temperature'][120:2640, :]
            zero = np.squeeze(zero)
            zero = reorder(zero)
            zero98[0, :] += np.mean(zero[0:3, :], axis=0)
            zero98[1, :] += np.mean(zero[3:6, :], axis=0)
            zero98[2, :] += np.mean(zero[6:9, :], axis=0)
            zero98[3, :] += np.mean(zero[9:12, :], axis=0)
            print('Done ' + d)
    zero98 = zero98/10
    for k in range(4):
        zero98[k, :] = np.ma.masked_array(zero98[k, :],
                                          np.logical_not(lsm[144:0:-1, :]))

    zero98b = np.zeros((4, 144, 192))
    for d in os.listdir('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/'):
        if d.find('0p98_') == 0:
            nc_f = d
            nc_fid = Dataset('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/' + nc_f, 'r')
            zerob = nc_fid.variables['air_temperature'][120:2640, :]
            zerob = np.squeeze(zerob)
            zerob = reorder(zerob)
            zero98b[0, :] += np.mean(zerob[0:3, :], axis=0)
            zero98b[1, :] += np.mean(zerob[3:6, :], axis=0)
            zero98b[2, :] += np.mean(zerob[6:9, :], axis=0)
            zero98b[3, :] += np.mean(zerob[9:12, :], axis=0)
            print('Done ' + d)
    zero98b = zero98b/10
    for k in range(4):
        zero98b[k, :] = np.ma.masked_array(zero98b[k, :],
                                           np.logical_not(lsm[144:0:-1, :]))

    one48 = np.zeros((4, 144, 192))
    for d in os.listdir('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/'):
        if d.find('1p48.') == 0:
            nc_f = d
            nc_fid = Dataset('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/' + nc_f, 'r')
            one = nc_fid.variables['air_temperature'][120:2640, :]
            one = np.squeeze(one)
            one = reorder(one)
            one48[0, :] += np.mean(one[0:3, :], axis=0)
            one48[1, :] += np.mean(one[3:6, :], axis=0)
            one48[2, :] += np.mean(one[6:9, :], axis=0)
            one48[3, :] += np.mean(one[9:12, :], axis=0)
            print('Done ' + d)
    one48 = one48/8
    for k in range(4):
        one48[k, :] = np.ma.masked_array(one48[k, :],
                                         np.logical_not(lsm[144:0:-1, :]))

    one48b = np.zeros((4, 144, 192))
    for d in os.listdir('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/'):
        if d.find('1p48_') == 0:
            nc_f = d
            nc_fid = Dataset('/home/bakerh/Documents/DPhil/CPDN/HAPPI/Klingaman/' + nc_f, 'r')
            oneb = nc_fid.variables['air_temperature'][120:2640, :]
            oneb = np.squeeze(oneb)
            oneb = reorder(oneb)
            one48b[0, :] += np.mean(oneb[0:3, :], axis=0)
            one48b[1, :] += np.mean(oneb[3:6, :], axis=0)
            one48b[2, :] += np.mean(oneb[6:9, :], axis=0)
            one48b[3, :] += np.mean(oneb[9:12, :], axis=0)
            print('Done ' + d)
    one48b = one48b/8
    for k in range(4):
        one48b[k, :] = np.ma.masked_array(one48b[k, :],
                                          np.logical_not(lsm[144:0:-1, :]))

    return control, zero98, zero98b, one48, one48b


def mapplotanom(plotdata, lat, lon, title=''):
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
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')
    lsm = lsm[0, 0, :]
    pdata = np.ma.masked_array(plotdata, mask=np.logical_not(lsm[144:0:-1, :]))
    pdata, lon = shiftgrid(180., pdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)

    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    my_cmap = plt.cm.jet(np.arange(256))
    my_cmap[110:146, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(-3, 3, 13)
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()
