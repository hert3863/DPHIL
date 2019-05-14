# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 13:50:41 2017

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


def processallmonthly(field):
    import glob
    from netCDF4 import Dataset
    output = {}
    exps = ['present', 'onepoint5',
            'onepoint5upper', 'onepoint5lower']
    for exp in exps:
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + field + '/*')
        data = np.zeros((4, 145, 192))
        for d in a:
            nc_fid = Dataset(d, 'r')
            na = nc_fid.variables[field][:, 0, :]
            data[0, :] += np.mean(np.concatenate((na[0:2, :],
                                  np.expand_dims(na[-1, :], axis=0)),
                                                 axis=0), axis=0)
            for j in range(10):
                data[0, :] += np.mean(na[12*j+11:12*j+14, :], axis=0)
            for i in range(11):
                data[1, :] += np.mean(na[12*i+2:12*i+5, :], axis=0)
                data[2, :] += np.mean(na[12*i+5:12*i+8, :], axis=0)
                data[3, :] += np.mean(na[12*i+8:12*i+11, :], axis=0)
            print('Done ' + d)
        data = data/(11*np.ma.size(a))
        output[exp] = data
    return output


def processalldaily(field):
    import glob
    from netCDF4 import Dataset
    output = {}
    exps = ['present', 'onepoint5',
            'onepoint5upper', 'onepoint5lower']
    for exp in exps:
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + field + '/*')
        data = np.zeros((4, 145, 192))
        for d in a:
            nc_fid = Dataset(d, 'r')
            na = nc_fid.variables[field][:, 0, :]
            data[0, :] += np.mean(np.concatenate((na[0:60, :],
                                  na[-30:, :]), axis=0), axis=0)
            for j in range(10):
                data[0, :] += np.mean(na[360*j+330:360*j+420, :], axis=0)
            for i in range(11):
                data[1, :] += np.mean(na[360*i+60:360*i+150, :], axis=0)
                data[2, :] += np.mean(na[360*i+150:360*i+240, :], axis=0)
                data[3, :] += np.mean(na[360*i+240:360*i+330, :], axis=0)
            print('Done ' + d)
        data = data/(11*np.ma.size(a))
        output[exp] = data
    return output


def processmix(itemcode):
    import glob
    from netCDF4 import Dataset
    output = {}
    exps = ['natall', 'natsst', 'natghg', 'present', 'onepoint5',
            'onepoint5upper', 'onepoint5lower']
    for exp in exps:
        if exp == 'present' or exp == 'onepoint5':
            a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                          '/atmos/' + itemcode + '_daily_mean/*')
            data = np.zeros((4, 144, 192))
            for d in a:
                nc_fid = Dataset(d, 'r')
                na = nc_fid.variables[itemcode + '_daily_mean'][:, 0, :]
                data[0, :] += np.mean(np.concatenate((na[0:60, :],
                                      na[-30:, :]), axis=0), axis=0)
                for j in range(10):
                    data[0, :] += np.mean(na[360*j+330:360*j+420, :], axis=0)
                for i in range(11):
                    data[1, :] += np.mean(na[360*i+60:360*i+150, :], axis=0)
                    data[2, :] += np.mean(na[360*i+150:360*i+240, :], axis=0)
                    data[3, :] += np.mean(na[360*i+240:360*i+330, :], axis=0)
                print('Done ' + d)
            data = data/(11*np.ma.size(a))
            output[exp] = data
        else:
            a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                          '/atmos/' + itemcode + '_monthly_mean/*')
            data = np.zeros((4, 144, 192))
            for d in a:
                nc_fid = Dataset(d, 'r')
                na = nc_fid.variables[itemcode + '_monthly_mean'][:, 0, :]
                data[0, :] += np.mean(np.concatenate((na[0:2, :],
                                      np.expand_dims(na[-1, :], axis=0)),
                                                     axis=0), axis=0)
                for j in range(10):
                    data[0, :] += np.mean(na[12*j+11:12*j+14, :], axis=0)
                for i in range(11):
                    data[1, :] += np.mean(na[12*i+2:12*i+5, :], axis=0)
                    data[2, :] += np.mean(na[12*i+5:12*i+8, :], axis=0)
                    data[3, :] += np.mean(na[12*i+8:12*i+11, :], axis=0)
                print('Done ' + d)
            data = data/(11*np.ma.size(a))
            output[exp] = data
    return output


def mapplotanom(pdata, lat, lon, mask='no', title=''):
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
    plt.figure()
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        pdata = np.ma.masked_array(pdata, mask=np.logical_not(lsm))
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
    ctrs = np.linspace(-1, 1, 13)
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    parallels = m.drawparallels(np.arange(-90., 91., 10.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def colourscaleint(m):
    """
    Takes data being plotted and normalises the colourscale between largest
    data value and it's negative multiple

    Parameters
    ----------
    m: float
        max of colourbar

    Returns
    -------
    caxismax: float
        max magnitude data value
    caxismin: float
        negative of max mag data value
    ctrs: array
        gradated colour scale contour array
    """
    m = abs(m)
    ctrs1 = np.arange(-m, 0, 0.05*m)
    ctrs2 = np.arange(0.05*m, 1.05*m, 0.05*m)
    ctrs = np.concatenate((ctrs1, ctrs2))
    caxismin = -m
    caxismax = m
    return caxismin, caxismax, ctrs


def mapplotprecip(pdata, lat, lon, mask='no', title=''):
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
    plt.figure()
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        pdata = np.ma.masked_array(pdata, mask=np.logical_not(lsm))
    pdata, lon = shiftgrid(180., pdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)

    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    mycmap2 = plt.cm.YlOrRd_r(np.arange(256))
    mycmap1 = plt.cm.Blues(np.arange(256))
    my_cmap = np.concatenate((mycmap2, mycmap1), axis=0)
    #my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = 0,.4,np.linspace(0,.4,21)
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def Tdist():
    field = 'item3236_daily_maximum'
    import glob
    from netCDF4 import Dataset
    output = {}
    exps = ['natall', 'natsst', 'natghg', 'present', 'onepoint5',
            'onepoint5upper', 'onepoint5lower']
    for exp in exps:
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + field + '/*')
        data = np.zeros((90*11*np.ma.size(a), 17, 28))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            na = nc_fid.variables[field][:, 0, 24:41, :28]
            for i in range(11):
                data[990*e+90*i:990*e+90*i+90] = na[360*i+150:360*i+240, :]
            print('Done ' + d)
        output[exp] = data
    return output


def Tdistmonth():
    field = 'item3236_monthly_mean'
    import glob
    from netCDF4 import Dataset
    output = {}
    exps = ['present', 'onepoint5',
            'onepoint5upper', 'onepoint5lower']
    for exp in exps:
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + field + '/*')
        data = np.zeros((3*11*np.ma.size(a), 17, 28))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            na = nc_fid.variables[field][:, 0, 24:41, :28]
            for i in range(11):
                data[33*e+3*i:33*e+3*i+3] = na[12*i+5:12*i+8, :]
            print('Done ' + d)
        output[exp] = data
    return output


def importallmonthly(field):
    import glob
    from netCDF4 import Dataset
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    masker=np.zeros((145,192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                masker[i, j] = 1
    output = {}
    exps = ['present', 'onepoint5', 'onepoint5upper', 'onepoint5lower']
    for exp in exps:
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + field + '/*')
        data = np.zeros((3*11*np.ma.size(a), 9391))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            na = nc_fid.variables[field][:, 0, :]
            na = na[:, masker == 0]
            for i in range(3):
                data[33*e+11*i:33*e+11*i+11, :] = na[5+i::12, :]
            print('Done ' + d)
        data = data[data > 0] - 273.15
        ey_data, ex_data = calc_return_times(data, 'descending')
        #conf_int = calc_return_time_confidences(data, direction='descending',
                                                #bsn=1000)
        output[exp] = ey_data, ex_data#, conf_int
    return output


def seanan(data, region):
    # region is array north to south extent, east to west extent
    data1 = np.copy(data)[:, region[0]:region[1], region[2]:region[3]]
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, region[0]:region[1], region[2]:region[3]]
    for i in range(np.ma.size(data, axis=0)):
        for j in range(region[1]-region[0]):
            for k in range(region[3]-region[2]):
                    if lsm[j, k] == 0:
                        data1[i, j, k] = np.nan
    return data1


def dailyTmaxreturn(field, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    output = {}
    exps = ['present', 'onepoint5', 'onepoint5upper', 'onepoint5lower']
    for exp in exps:
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + field + '/*')
        data = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            na = nc_fid.variables[field][:, 0, region[0]:region[1],
                                         region[2]:region[3]]
            for i in range(11):
                data[990*e+90*i:990*e+90*i+90, :] = na[360*i+150:360*i+240, :]
            print('Done ' + d)
        data = data[data > 0] - 273.15
        ey_data, ex_data = calc_return_times(data, 'descending')
        #conf_int = calc_return_time_confidences(data, direction='descending',
                                                #bsn=1000)
        output[exp] = ey_data, ex_data#, conf_int
    return output
