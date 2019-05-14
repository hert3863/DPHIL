# -*- coding: utf-8 -*-
"""
Created on Tue Apr 25 15:05:16 2017

@author: bakerh
"""

import numpy as np


def plotall(data, colorlimit, meaning='mean', mask='yes', precip='no'):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_522/atmos/item15201_monthly_mean/item15201_monthly_mean_a00a_2090-01_2100-12.nc', 'latitude1')
    if meaning == 'mean':
        winter = np.mean(data['batch_522'][0, :], axis=0) - np.mean(data['batch_521'][0, :], axis=0)
        summer = np.mean(data['batch_522'][1, :], axis=0) - np.mean(data['batch_521'][1, :], axis=0)
        sig0 = ttest(data['batch_522'][0, :],
                     data['batch_521'][0, :], siglevel=10)
        sig1 = ttest(data['batch_522'][1, :],
                     data['batch_521'][1, :], siglevel=10)
    elif meaning == 'non':
        winter = data['batch_522'][0, :] - data['batch_521'][0, :]
        summer = data['batch_522'][1, :] - data['batch_521'][1, :]
        sig0 = np.zeros((145, 192))
        sig1 = np.zeros((145, 192))
    else:
        winter = np.mean(data['batch_522'][meaning][0, :], axis=0) - np.mean(data['batch_521'][meaning][0, :], axis=0)
        summer = np.mean(data['batch_522'][meaning][1, :], axis=0) - np.mean(data['batch_521'][meaning][1, :], axis=0)
        sig0 = ttest(data['batch_522'][meaning][0, :],
                     data['batch_521'][meaning][0, :], siglevel=10)
        sig1 = ttest(data['batch_522'][meaning][1, :],
                     data['batch_521'][meaning][1, :], siglevel=10)
    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        winter = np.ma.masked_array(winter, mask=np.logical_not(lsm))
        summer = np.ma.masked_array(summer, mask=np.logical_not(lsm))
        sig0 = np.ma.masked_array(sig0, mask=np.logical_not(lsm))
        sig1 = np.ma.masked_array(sig1, mask=np.logical_not(lsm))

    winter, lon1 = shiftgrid(180., winter, lon, start=False)
    winter, lon1 = addcyclic(winter, lon1)
    sig0, lon2 = shiftgrid(180., sig0, lon, start=False)
    sig0, lon2 = addcyclic(sig0, lon2)
    summer, lon3 = shiftgrid(180., summer, lon, start=False)
    summer, lon3 = addcyclic(summer, lon3)
    sig1, lon4 = shiftgrid(180., sig1, lon, start=False)
    sig1, lon4 = addcyclic(sig1, lon4)
    meshlon, meshlat = np.meshgrid(lon1, lat)
    ctrs = np.linspace(-colorlimit, colorlimit, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    if precip == 'yes':
        my_cmap = my_cmap[::-1]
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap

    ax1 = fig.add_subplot(1, 2, 1)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, winter, ctrs,
                      cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    plt.contourf(x, y, sig0, levels=[0, .999999, 1.9, 2.1],
                 hatches=["", "////", "xx"], alpha=0)
    ax2 = fig.add_subplot(1, 2, 2)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot1 = m.contourf(x, y, summer, ctrs,
                       cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                       extend='both')
    # parallels = m.drawparallels(np.arange(-90., 91., 15.))
    # meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    # m.drawparallels(parallels, labels=[True, True, True, True])
    # m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.contourf(x, y, sig1, levels=[0, .999999, 1.9, 2.1],
                 hatches=["", "////", "xx"], alpha=0)

    cbar_ax = fig.add_axes([0.2, 0.3, 0.6, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max')
    b.set_label(label='Change in WBGT_95% [K]', size=20,
                weight='bold', fontsize=20)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)

    pad = 15  # in points
    ax1.annotate('DJF', xy=(0.5, 1), xytext=(0, pad), fontweight='bold',
                 fontsize=32,
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('JJA', xy=(0.5, 1), xytext=(0, pad), fontweight='bold',
                 fontsize=32,
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax1.annotate('WBGT_95%', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontweight='bold', fontsize=32,
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    plt.show()
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.15, left=.05,
                        right=.95)


def plotall2(data, data2, colorlimit, meaning='mean', meaning2='mean',
             mask='yes', precip='no'):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_522/atmos/item15201_monthly_mean/item15201_monthly_mean_a00a_2090-01_2100-12.nc', 'latitude1')
    if meaning == 'mean':
        winter = np.mean(data['batch_522'][0, :], axis=0) - np.mean(data['batch_521'][0, :], axis=0)
        summer = np.mean(data['batch_522'][1, :], axis=0) - np.mean(data['batch_521'][1, :], axis=0)
        sig0 = ttest(data['batch_522'][0, :],
                     data['batch_521'][0, :], siglevel=10)
        sig1 = ttest(data['batch_522'][1, :],
                     data['batch_521'][1, :], siglevel=10)
    else:
        winter = np.mean(data['batch_522'][meaning][0, :], axis=0) - np.mean(data['batch_521'][meaning][0, :], axis=0)
        summer = np.mean(data['batch_522'][meaning][1, :], axis=0) - np.mean(data['batch_521'][meaning][1, :], axis=0)
        sig0 = ttest(data['batch_522'][meaning][0, :],
                     data['batch_521'][meaning][0, :], siglevel=10)
        sig1 = ttest(data['batch_522'][meaning][1, :],
                     data['batch_521'][meaning][1, :], siglevel=10)
    if meaning2 == 'mean':
        winter2 = np.mean(data2['batch_522'][0, :], axis=0) - np.mean(data2['batch_521'][0, :], axis=0)
        summer2 = np.mean(data2['batch_522'][1, :], axis=0) - np.mean(data2['batch_521'][1, :], axis=0)
        sig3 = ttest(data2['batch_522'][0, :],
                     data2['batch_521'][0, :], siglevel=10)
        sig4 = ttest(data2['batch_522'][1, :],
                     data2['batch_521'][1, :], siglevel=10)
    else:
        winter2 = np.mean(data2['batch_522'][meaning2][0, :], axis=0) - np.mean(data2['batch_521'][meaning2][0, :], axis=0)
        summer2 = np.mean(data2['batch_522'][meaning2][1, :], axis=0) - np.mean(data2['batch_521'][meaning2][1, :], axis=0)
        sig3 = ttest(data2['batch_522'][meaning2][0, :],
                     data2['batch_521'][meaning2][0, :], siglevel=10)
        sig4 = ttest(data2['batch_522'][meaning2][1, :],
                     data2['batch_521'][meaning2][1, :], siglevel=10)
    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        winter = np.ma.masked_array(winter, mask=np.logical_not(lsm))
        summer = np.ma.masked_array(summer, mask=np.logical_not(lsm))
        sig0 = np.ma.masked_array(sig0, mask=np.logical_not(lsm))
        sig1 = np.ma.masked_array(sig1, mask=np.logical_not(lsm))
        winter2 = np.ma.masked_array(winter2, mask=np.logical_not(lsm))
        summer2 = np.ma.masked_array(summer2, mask=np.logical_not(lsm))
        sig3 = np.ma.masked_array(sig3, mask=np.logical_not(lsm))
        sig4 = np.ma.masked_array(sig4, mask=np.logical_not(lsm))

    winter, lon1 = shiftgrid(180., winter, lon, start=False)
    winter, lon1 = addcyclic(winter, lon1)
    sig0, lon2 = shiftgrid(180., sig0, lon, start=False)
    sig0, lon2 = addcyclic(sig0, lon2)
    summer, lon3 = shiftgrid(180., summer, lon, start=False)
    summer, lon3 = addcyclic(summer, lon3)
    sig1, lon4 = shiftgrid(180., sig1, lon, start=False)
    sig1, lon4 = addcyclic(sig1, lon4)
    winter2, lon5 = shiftgrid(180., winter2, lon, start=False)
    winter2, lon5 = addcyclic(winter2, lon5)
    sig3, lon6 = shiftgrid(180., sig3, lon, start=False)
    sig3, lon6 = addcyclic(sig3, lon6)
    summer2, lon7 = shiftgrid(180., summer2, lon, start=False)
    summer2, lon7 = addcyclic(summer2, lon7)
    sig4, lon8 = shiftgrid(180., sig4, lon, start=False)
    sig4, lon8 = addcyclic(sig4, lon8)
    meshlon, meshlat = np.meshgrid(lon1, lat)
    ctrs = np.linspace(-colorlimit, colorlimit, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    if precip == 'yes':
        my_cmap = my_cmap[::-1]
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap

    ax1 = fig.add_subplot(2, 2, 1)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, winter, ctrs,
                      cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    plt.contourf(x, y, sig0, levels=[0, .999999, 1.9, 2.1],
                 hatches=["", "////", "xx"], alpha=0)
    ax2 = fig.add_subplot(2, 2, 2)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, summer, ctrs,
               cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
               extend='both')
    plt.contourf(x, y, sig1, levels=[0, .999999, 1.9, 2.1],
                 hatches=["", "////", "xx"], alpha=0)
    ax3 = fig.add_subplot(2, 2, 3)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, winter2, ctrs,
               cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
               extend='both')
    plt.contourf(x, y, sig3, levels=[0, .999999, 1.9, 2.1],
                 hatches=["", "////", "xx"], alpha=0)
    ax4 = fig.add_subplot(2, 2, 4)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, summer2, ctrs,
               cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs), extend='both')
    plt.contourf(x, y, sig4, levels=[0, .999999, 1.9, 2.1],
                 hatches=["", "////", "xx"], alpha=0)
    cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max')
    b.set_label(label='Temperature difference [K]', size=20,
                weight='bold', fontsize=20)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)

    pad = 15  # in points
    ax1.annotate('DJF', xy=(0.5, 1), xytext=(0, pad), fontweight='bold',
                 fontsize=32,
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('JJA', xy=(0.5, 1), xytext=(0, pad), fontweight='bold',
                 fontsize=32,
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax1.annotate('TX', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontweight='bold', fontsize=32,
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('TXx', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontweight='bold', fontsize=32,
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    plt.show()
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.15, left=.05,
                        right=.95)


def ttest(series1, series2, siglevel=10, testtype='two'):
    """
    Student t-test for a time series

    Parameters
    ----------
    series1: array
        control run
    series2: array
        forced run

    Returns
    -------
    sig: bool
        is it significant
    sig: bool
        is it discernible and significant
    """
    import scipy.stats as st
    sigarray = np.full((145, 192), siglevel)
    sig = np.zeros((145, 192))
    if testtype == 'two':
        a = 2
    elif testtype == 'one':
        a = 1
    else:
        print("Error, test type must be 'one' or 'two'")
    z, p = st.ttest_ind(series1, series2, axis=0, equal_var=False)
    d = np.sqrt(np.var(series1, axis=0, ddof=1) + np.var(series2, axis=0,
                                                         ddof=1))
    z1 = (np.mean(series1, axis=0) - np.mean(series2, axis=0)) / d
    p1 = 1 - st.norm.cdf(np.abs(z1))
    sig = np.greater_equal(sigarray, p*100*a).astype(int)
    sig_d = np.greater_equal(sigarray, p1*100*a).astype(int)
    sig_art = sig + sig_d
    return sig_art


def monthmean(item, level, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*3))
    windices = np.zeros((11*3))
    for i in range(11):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    for i in range(11):
        windices[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}
    output_95 = {}
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + item + '/*')
        tbar = np.zeros((2, np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        sdata = np.zeros((3*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((3*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        t_95 = np.zeros((2, 145, 192))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[33*e:33*(e+1)] = nc_fid.variables[item
                                                    ][sindices, level,
                                                      region[0]:
                                                      region[1],
                                                      region[2]:region[3]]
            wdata[33*e:33*(e+1)] = nc_fid.variables[item
                                                    ][windices, level,
                                                      region[0]:
                                                      region[1],
                                                      region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            tbar[0, y, :] = np.mean(wdata[33*y:33*(y+1), :], axis=0)
            tbar[1, y, :] = np.mean(sdata[33*y:33*(y+1), :], axis=0)
        wdataTmonth = np.zeros((np.ma.size(a)*11, 145, 192))
        sdataTmonth = np.zeros((np.ma.size(a)*11, 145, 192))
        for z in range(np.ma.size(a)):
            wdataTmonth[z, :] = np.mean(wdata[3*z:3*(z+1), :], axis=0)
            sdataTmonth[z, :] = np.mean(sdata[3*z:3*(z+1), :], axis=0)
        valuew = np.percentile(wdataTmonth, 95, axis=0)
        values = np.percentile(sdataTmonth, 95, axis=0)
        for i in range(145):
            for j in range(192):
                t_gridw = wdataTmonth[:, i, j]
                t_95[0, i, j] = np.mean(t_gridw[wdataTmonth[:, i, j] >
                                        valuew[i, j]])
                t_grids = sdataTmonth[:, i, j]
                t_95[1, i, j] = np.mean(t_grids[sdataTmonth[:, i, j] >
                                        values[i, j]])
        output[exp] = tbar
        output_95[exp] = t_95
    return output, output_95


def monthmeanall(item, level, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']
    means = np.zeros((4))
    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + item + '/*')
        b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + 'item2205_monthly_mean' + '/*')
        data = np.zeros((132*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        data1 = np.zeros((132*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            nc_fib = Dataset(b[e], 'r')
            data[132*e:132*(e+1)] = nc_fid.variables[item
                                                     ][:, level,
                                                       region[0]:
                                                       region[1],
                                                       region[2]:region[3]]
            data1[132*e:132*(e+1)] = nc_fib.variables['item2205_monthly_mean'
                                                      ][:, level,
                                                        region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        data = data + data1
        data = np.mean(data, axis=0)
        means[x] = latweightmean(data)

    return means


def dailymean(item, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}

    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/' + item + '/*')
        tbar = np.zeros((2, np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables[item
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[990*e:990*(e+1)] = nc_fid.variables[item
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            tbar[0, y, :] = np.mean(wdata[990*y:990*(y+1), :], axis=0)
            tbar[1, y, :] = np.mean(sdata[990*y:990*(y+1), :], axis=0)

        output[exp] = tbar
    return output


def precipall(region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/batch_518' +
                  '/atmos/item5216_daily_mean/*')
    sdatar90p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    wdatar90p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatar90p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        wdatar90p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        print('Done R90p calc ' + str(e+1))
    days = np.arange(0, 990*np.ma.size(a), 90).astype(int)
    r90 = np.zeros((2, 90, region[1]-region[0],
                    region[3]-region[2]))
    for i in range(90):
        r90[0, i, :] = np.percentile(wdatar90p[days+i, :], 90, axis=0)
        r90[1, i, :] = np.percentile(sdatar90p[days+i, :], 90, axis=0)

    outputp = {}
    outputrx1d = {}
    outputr90p = {}
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item5216_daily_mean/*')
        pbar = np.zeros((2, np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        rx1d = np.zeros((2, np.ma.size(a), 11, region[1]-region[0],
                        region[3]-region[2]))
        r90p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                        region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done precip ' + str(exp) + ' ' + str(e+1))
        days = np.arange(0, 990, 90).astype(int)
        truth = np.zeros((2, 90*11*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        for k in range(11*np.ma.size(a)):
            truth[0, 90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                          r90[0, :]).astype(int)
            truth[1, 90*k:90*(k+1), :] = (sdata[90*k:90*(k+1), :] >
                                          r90[1, :]).astype(int)
        for m in range(np.ma.size(a)):
            r90p[0, m, :] = np.sum(truth[0, 990*m:990*(m+1), :], axis=0)
            r90p[1, m, :] = np.sum(truth[1, 990*m:990*(m+1), :], axis=0)
        r90p = r90p / (9.9)  # percentage
        for y in range(np.ma.size(a)):
            pbar[0, y, :] = np.mean(wdata[990*y:990*(y+1), :], axis=0)
            pbar[1, y, :] = np.mean(sdata[990*y:990*(y+1), :], axis=0)
            for t in range(11):
                rx1d[0, y, t, :] = np.max(wdata[990*y+90*t:990*y+90*(t+1), :],
                                          axis=0)
                rx1d[1, y, t, :] = np.max(sdata[990*y+90*t:990*y+90*(t+1), :],
                                          axis=0)
        rx1d = np.mean(rx1d, axis=2)
        outputp[exp] = pbar * 86400
        outputrx1d[exp] = rx1d * 86400
        outputr90p[exp] = r90p
    return outputp, outputrx1d, outputr90p


def R95ptot(region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/batch_518' +
                  '/atmos/item5216_daily_mean/*')
    sdatar95p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    wdatar95p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatar95p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        wdatar95p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        print('Done R95p calc ' + str(e+1))
    r95 = np.ones((2, region[1]-region[0], region[3]-region[2]))
    for i in range(region[1]-region[0]):
        for j in range(region[3]-region[2]):
            wdatar95p_d = wdatar95p[:, i, j] * 86400
            sdatar95p_d = sdatar95p[:, i, j] * 86400
            wdatar95p_d1 = wdatar95p_d[wdatar95p_d > 1]
            sdatar95p_d1 = sdatar95p_d[sdatar95p_d > 1]
            if np.ma.size(wdatar95p_d1) != 0:
                r95[0, i, j] = np.percentile(wdatar95p_d1, 95)
            if np.ma.size(sdatar95p_d1) != 0:
                r95[1, i, j] = np.percentile(sdatar95p_d1, 95)

    outputr95p = {}
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item5216_daily_mean/*')

        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        r95p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                        region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            wdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            print('Done precip ' + str(exp) + ' ' + str(e+1))
        for m in range(np.ma.size(a)):
            r95p[0, m] = np.sum(wdata[990*m:990*(m+1), :] >
                                r95[0], axis=0)
            r95p[1, m] = np.sum(sdata[990*m:990*(m+1), :] >
                                r95[1], axis=0)
        outputr95p[exp] = r95p / 11  # annual precip change
    return outputr95p


def tmax(region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/batch_518' +
                  '/atmos/item3236_daily_maximum/*')
    sdatatx90p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    wdatatx90p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatatx90p[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                       ][sindices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        wdatatx90p[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                       ][windices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        print('Done TX90p calc ' + str(e+1))
    days = np.arange(0, 990*np.ma.size(a), 90).astype(int)
    t90 = np.zeros((2, 90, region[1]-region[0],
                    region[3]-region[2]))
    for i in range(90):
        t90[0, i, :] = np.percentile(wdatatx90p[days+i, :], 90, axis=0)
        t90[1, i, :] = np.percentile(sdatatx90p[days+i, :], 90, axis=0)

    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']
    output = {}

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item3236_daily_maximum/*')
        tx90p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                          region[3]-region[2]))
        txx = np.zeros((2, np.ma.size(a), 11, region[1]-region[0],
                        region[3]-region[2]))
        tx = np.zeros((2, np.ma.size(a), region[1]-region[0],
                       region[3]-region[2]))
        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done txmax ' + str(exp) + ' ' + str(e+1))

        days = np.arange(0, 990, 90).astype(int)
        truth = np.zeros((2, 90*11*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        for k in range(11*np.ma.size(a)):
            truth[0, 90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                          t90[0, :]).astype(int)
            truth[1, 90*k:90*(k+1), :] = (sdata[90*k:90*(k+1), :] >
                                          t90[1, :]).astype(int)
        for m in range(np.ma.size(a)):
            tx90p[0, m, :] = np.sum(truth[0, 990*m:990*(m+1), :], axis=0)
            tx90p[1, m, :] = np.sum(truth[1, 990*m:990*(m+1), :], axis=0)
        tx90p = tx90p / (9.9)  # percentage

        for y in range(np.ma.size(a)):
            tx[0, y, :] = np.mean(wdata[990*y:990*(y+1), :], axis=0)
            tx[1, y, :] = np.mean(sdata[990*y:990*(y+1), :], axis=0)
            for t in range(11):
                txx[0, y, t, :] = np.max(wdata[990*y+90*t:990*y+90*(t+1), :],
                                         axis=0)
                txx[1, y, t, :] = np.max(sdata[990*y+90*t:990*y+90*(t+1), :],
                                         axis=0)
        txx = np.mean(txx, axis=2)
        output[exp] = tx, txx, tx90p
    return output


def tmin(region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/batch_518' +
                  '/atmos/item3236_daily_minimum/*')
    sdatatn10p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    wdatatn10p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatatn10p[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_minimum'
                                                       ][sindices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        wdatatn10p[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_minimum'
                                                       ][windices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        print('Done ' + str(e+1))
    days = np.arange(0, 990*np.ma.size(a), 90).astype(int)
    t10 = np.zeros((2, 90, region[1]-region[0],
                    region[3]-region[2]))
    for i in range(90):
        t10[0, i, :] = np.percentile(wdatatn10p[days+i, :], 10, axis=0)
        t10[1, i, :] = np.percentile(sdatatn10p[days+i, :], 10, axis=0)

    exps = ['batch_520', 'batch_521', 'batch_522']
    output = {}

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item3236_daily_minimum/*')
        tn10p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                          region[3]-region[2]))
        tnn = np.zeros((2, np.ma.size(a), 11, region[1]-region[0],
                        region[3]-region[2]))
        tn = np.zeros((2, np.ma.size(a), region[1]-region[0],
                       region[3]-region[2]))
        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_minimum'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_minimum'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done ' + str(exp) + ' ' + str(e+1))

        days = np.arange(0, 990, 90).astype(int)
        truth = np.zeros((2, 90*11*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        for k in range(11*np.ma.size(a)):
            truth[0, 90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] <
                                          t10[0, :]).astype(int)
            truth[1, 90*k:90*(k+1), :] = (sdata[90*k:90*(k+1), :] <
                                          t10[1, :]).astype(int)
        for m in range(np.ma.size(a)):
            tn10p[0, m, :] = np.sum(truth[0, 990*m:990*(m+1), :], axis=0)
            tn10p[1, m, :] = np.sum(truth[1, 990*m:990*(m+1), :], axis=0)
        tn10p = tn10p / (9.9)  # percentage

        for y in range(np.ma.size(a)):
            tn[0, y, :] = np.mean(wdata[990*y:990*(y+1), :], axis=0)
            tn[1, y, :] = np.mean(sdata[990*y:990*(y+1), :], axis=0)
            for t in range(11):
                tnn[0, y, t, :] = np.min(wdata[990*y+90*t:990*y+90*(t+1), :],
                                         axis=0)
                tnn[1, y, t, :] = np.min(sdata[990*y+90*t:990*y+90*(t+1), :],
                                         axis=0)
        tnn = np.mean(tnn, axis=2)
        output[exp] = tn, tnn, tn10p
    return output


def wbgt_calc():
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*3))
    windices = np.zeros((11*3))
    for i in range(11):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    for i in range(11):
        windices[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}
    output_95 = {}
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item3236_monthly_mean/*')
        wbgt = np.zeros((2, np.ma.size(a), 145, 192))
        wbgt_95 = np.zeros((2, 145, 192))
        sdataT = np.zeros((3*11*np.ma.size(a), 145, 192))
        wdataT = np.zeros((3*11*np.ma.size(a), 145, 192))
        sdataH = np.zeros((3*11*np.ma.size(a), 145, 192))
        wdataH = np.zeros((3*11*np.ma.size(a), 145, 192))
        wbgts = np.zeros((3*11*np.ma.size(a), 145, 192))
        wbgtw = np.zeros((3*11*np.ma.size(a), 145, 192))
        for e, d in enumerate(a):
            b = d[:58]+str(45)+d[60:80]+str(45)+d[82:]
            nc_fid = Dataset(d, 'r')
            nc_fib = Dataset(b, 'r')

            sdataT[33*e:33*(e+1)] = nc_fid.variables['item3236_monthly_mean'
                                                     ][sindices, 0, :]
            wdataT[33*e:33*(e+1)] = nc_fid.variables['item3236_monthly_mean'
                                                     ][windices, 0, :]
            sdataH[33*e:33*(e+1)] = nc_fib.variables['item3245_monthly_mean'
                                                     ][sindices, 0, :]
            wdataH[33*e:33*(e+1)] = nc_fib.variables['item3245_monthly_mean'
                                                     ][windices, 0, :]
            print('Done ' + str(exp) + ' ' + str(e+1))
        wdataTmonth = np.zeros((np.ma.size(a)*11, 145, 192))
        sdataTmonth = np.zeros((np.ma.size(a)*11, 145, 192))
        wdataHmonth = np.zeros((np.ma.size(a)*11, 145, 192))
        sdataHmonth = np.zeros((np.ma.size(a)*11, 145, 192))
        wwbgtmonth = np.zeros((np.ma.size(a)*11, 145, 192))
        swbgtmonth = np.zeros((np.ma.size(a)*11, 145, 192))

        Es = 6.112*np.exp(17.62*(sdataT-273.15)/(243.12+sdataT-273.15))
        Ew = 6.112*np.exp(17.62*(wdataT-273.15)/(243.12+wdataT-273.15))
        es = sdataH*Es/100
        ew = wdataH*Ew/100
        wbgts = 0.567*(sdataT-273.15)+0.393*es+3.94
        wbgtw = 0.567*(wdataT-273.15)+0.393*ew+3.94
        for y in range(np.ma.size(a)):
            wbgt[0, y, :] = np.mean(wbgtw[33*y:33*(y+1), :], axis=0)
            wbgt[1, y, :] = np.mean(wbgts[33*y:33*(y+1), :], axis=0)
        for z in range(np.ma.size(a)):
            wdataTmonth[z, :] = np.mean(wdataT[3*z:3*(z+1), :], axis=0)
            sdataTmonth[z, :] = np.mean(sdataT[3*z:3*(z+1), :], axis=0)
            wdataHmonth[z, :] = np.mean(wdataH[3*z:3*(z+1), :], axis=0)
            sdataHmonth[z, :] = np.mean(sdataH[3*z:3*(z+1), :], axis=0)
            wwbgtmonth[z, :] = np.mean(wbgtw[3*z:3*(z+1), :], axis=0)
            swbgtmonth[z, :] = np.mean(wbgts[3*z:3*(z+1), :], axis=0)
        valuew = np.percentile(wdataTmonth, 95, axis=0)
        values = np.percentile(sdataTmonth, 95, axis=0)
        for i in range(145):
            for j in range(192):
                wbgt_gridw = wwbgtmonth[:, i, j]
                wbgt_95[0, i, j] = np.mean(wbgt_gridw[wdataTmonth[:, i, j] >
                                           valuew[i, j]])
                wbgt_grids = swbgtmonth[:, i, j]
                wbgt_95[1, i, j] = np.mean(wbgt_grids[sdataTmonth[:, i, j] >
                                           values[i, j]])
        output[exp] = wbgt
        output_95[exp] = wbgt_95
    return output, output_95


def latweightmean(data, msk='no'):
    from netcdfread import ncread
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    weighted = data * meshlatweight
    if msk == 'yes':
        mask = np.zeros((np.ma.size(data, axis=0), 145, 192))
        mask[:] = lsm1
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'NH':
        mask = np.zeros((np.ma.size(data, axis=0), 73, 192))
        lsm1 = lsm1[:73, :]
        mask[:] = lsm1
        weighted = weighted[:, :73, :]
        meshlatweight = meshlatweight[:73, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'TROP':
        mask = np.zeros((np.ma.size(data, axis=0), 49, 192))
        lsm1 = lsm1[48:97, :]
        mask[:] = lsm1
        weighted = weighted[:, 48:97, :]
        meshlatweight = meshlatweight[48:97, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'MID':
        mask = np.zeros((np.ma.size(data, axis=0), 48, 192))
        lsm1 = lsm1[:48, :]
        mask[:] = lsm1
        weighted = weighted[:,  :48, :]
        meshlatweight = meshlatweight[:48, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    else:
        meaned = np.mean(weighted) / np.mean(meshlatweight)
    return meaned


def totals():
    import glob
    from netCDF4 import Dataset
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']
    total_r = np.zeros((4, 2, 3))  # exp; global, land only;
    total_t = np.zeros((4, 2, 3))  # all, winter, summer
    percentiles = np.zeros((4, 3))
    globprcnts = np.zeros((4, 3))
    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item5216_daily_mean/*')
        b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item3236_monthly_mean/*')
        dataR = np.zeros((360*11*np.ma.size(a), 145, 192))
        dataT = np.zeros((12*11*np.ma.size(b), 145, 192))

        for e, d in enumerate(a):
            sindices = np.zeros((11*90*np.ma.size(a)))
            windices = np.zeros((11*90*np.ma.size(a)))
            for i in range(11*np.ma.size(a)):
                sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
            for i in range(11*np.ma.size(a)):
                windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i,
                                                          59+360*i, 60),
                                                         np.linspace(330+360*i,
                                                         359+360*i, 30)))
            sindices = sindices.astype(int)
            windices = windices.astype(int)

            nc_fid = Dataset(d, 'r')
            dataR[3960*e:3960*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                        ][:, 0, :]
            print('Done R totals ' + str(exp) + ' ' + str(e+1))
        total_r[x, 0, 0] = latweightmean(dataR) * 86400
        total_r[x, 0, 1] = latweightmean(dataR[windices, :]) * 86400
        total_r[x, 0, 2] = latweightmean(dataR[sindices, :]) * 86400
        total_r[x, 1, 0] = latweightmean(dataR, msk='yes') * 86400
        dataRw = dataR[windices, :]
        dataRs = dataR[sindices, :]
        total_r[x, 1, 1] = latweightmean(dataRw, msk='yes') * 86400
        total_r[x, 1, 2] = latweightmean(dataRs, msk='yes') * 86400

        for f, g in enumerate(b):
            sindicesm = np.zeros((11*3*np.ma.size(b)))
            windicesm = np.zeros((11*3*np.ma.size(b)))
            for i in range(11*np.ma.size(b)):
                sindicesm[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
            for i in range(11*np.ma.size(b)):
                windicesm[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
            sindicesm = sindicesm.astype(int)
            windicesm = windicesm.astype(int)
            nc_fid = Dataset(g, 'r')
            dataT[132*f:132*(f+1)] = nc_fid.variables['item3236_monthly_mean'
                                                      ][:, 0, :]
            print('Done T totals ' + str(exp) + ' ' + str(f+1))
        total_t[x, 0, 0] = latweightmean(dataT)
        total_t[x, 0, 1] = latweightmean(dataT[windicesm, :])
        total_t[x, 0, 2] = latweightmean(dataT[sindicesm, :])
        total_t[x, 1, 0] = latweightmean(dataT, msk='yes')
        dataTw = dataT[windicesm, :]
        dataTs = dataT[sindicesm, :]
        total_t[x, 1, 1] = latweightmean(dataTw, msk='yes')
        total_t[x, 1, 2] = latweightmean(dataTs, msk='yes')
        dataT_point = np.zeros((np.ma.size(b)*132))
        for l in range(np.ma.size(b)*132):
            dataT_point[l] = latweightmean(dataT[l, :])
        dataT_11 = np.zeros((np.ma.size(b)))
        for y in range(np.ma.size(b)):
            dataT_11[y] = np.mean(dataT_point[132*y:132*(y+1)], axis=0)
        percentiles[x, :] = np.percentile(dataT_11, [10, 50, 90])
        for a in range(4):
            globprcnts[a, 0] = percentiles[a, 1]
            globprcnts[a, 1] = percentiles[a, 1] - percentiles[a, 0]
            globprcnts[a, 2] = percentiles[a, 2] - percentiles[a, 1]
    return total_r, total_t, globprcnts


def percents(tbar, txmax, r, r90p, sw, lw, wbgt):
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']
    percentiles = np.zeros((4, 11, 3))
    prcnts = np.zeros((4, 11, 3))
    all_mean = np.zeros((4, 11))

    for e, exp in enumerate(exps):
        tbar1 = tbar[exp][1, :, :, :]
        tx1 = txmax[exp][0][1, :, :, :]
        txx1 = txmax[exp][1][1, :, :, :]
        tx90p1 = txmax[exp][2][1, :, :, :]
        sw1 = sw[exp][1, :, :, :]
        lw1 = lw[exp][1, :, :, :]
        wbgt1 = wbgt[exp][1, :, :, :]
        r1 = r[exp][1, :, :, :]
        r90p1 = r90p[exp][1, :, :, :]
        varis = [tbar1, tx1, txx1, tx90p1, wbgt1, sw1, lw1]
        varis_r = [r1, r90p1]
        for i, var in enumerate(varis):
            means = []
            for j in range(np.ma.size(var, axis=0)):
                means.append(latweightmean(np.expand_dims(var[j, :], axis=0), msk='NH'))
            percentiles[e, i, :] = np.percentile(means, [10, 50, 90])
            all_mean[e, i] = np.mean(means)
        for i, var in enumerate(varis_r):
            means = []
            means_tropical = []
            for j in range(np.ma.size(var, axis=0)):
                means.append(latweightmean(np.expand_dims(var[j, :], axis=0), msk='MID'))
                means_tropical.append(latweightmean(np.expand_dims(var[j, :], axis=0), msk='TROP'))
            percentiles[e, 2*i+7, :] = np.percentile(means, [10, 50, 90])
            percentiles[e, 2*i+8, :] = np.percentile(means_tropical,
                                                     [10, 50, 90])
            all_mean[e, 2*i+7] = np.mean(means)
            all_mean[e, 2*i+8] = np.mean(means_tropical)

    for a in range(4):
        for b in range(11):
            prcnts[a, b, 0] = percentiles[a, b, 1]
            prcnts[a, b, 1] = percentiles[a, b, 1] - percentiles[a, b, 0]
            prcnts[a, b, 2] = percentiles[a, b, 2] - percentiles[a, b, 1]
    return prcnts, all_mean


def main():
    tbar = monthmean('item3236_monthly_mean', 0)
    txmax = tmax()
    r, rx1d, r90p = precipall()
    lw = dailymean('item2207_daily_mean')
    sw = dailymean('item1235_daily_mean')
    total_r, total_t, globprcnts = totals()
    prcnts, all_mean = percents(tbar, txmax, r, r90p, sw, lw)
    return tbar, txmax, r, rx1d, r90p, lw, sw, total_r, total_t, globprcnts, prcnts, all_mean


def mlr(tmean, var, norm='no'):
    import pandas as pd
    from statsmodels.formula.api import ols

    co2 = np.array([390.4, 486.6, 395.8, 550.0])
    f = 3.74 * np.log(co2/278) / np.log(2)  # 278 as preindust [CO2]
    d_f = f - f[0]
    d_t = tmean-tmean[0]
    d_var = var-var[0]
    # d_co2 = co2-co2[0]
    # logco2 = np.log(co2)
    # d_logco2 = logco2-logco2[0]
    if norm == 'yes':
        d_t /= d_t[1]
        d_var /= d_var[1]
    data = pd.DataFrame({'d_f': d_f, 'd_t': d_t, 'd_var': d_var})
    model = ols("d_var ~ d_t + d_f - 1", data).fit()
    # print(model.summary())
    # print(model._results.params)
    return model._results.params


def fig4(tmean, var, conversion):
    from matplotlib import pyplot as plt

    def gtc2ppm(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        return ["%.0f" % z for z in x2]

    def gtc2ppmfloat(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        return x2

    def ppm2gtc(x1):
        x2 = np.interp(x1, conversion[:, 1], conversion[:, 2])
        return x2

    def df2gtc(x1):
        dco2 = 390.4*np.exp(np.log(2)*x1/3.74)
        x2 = np.interp(dco2, conversion[:, 1], conversion[:, 2])
        return x2

    def dfdc(c):
        ppm1 = gtc2ppmfloat(c-10)
        ppm2 = gtc2ppmfloat(c+10)
        f1 = 3.74 * np.log(ppm1/278) / np.log(2)
        f2 = 3.74 * np.log(ppm2/278) / np.log(2)
        dfdc = (f2-f1)/20
        return dfdc

    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    co2 = np.array([390.4, 486.6, 395.8, 550.0])
    beta, alpha = mlr(tmean, var)
    co2gtc = ppm2gtc(co2)
    g2 = 0.001
    g1 = 0.0025
    ax1.fill_between(conversion[:, 2], g2*conversion[:, 2],
                     g1*conversion[:, 2], facecolor='gray')
    ax1.scatter(co2gtc[1:], tmean[1:]-tmean[0], color='b', marker='o', s=20)
    ax1.plot(np.arange(0, 1010, 10), -alpha*dfdc(np.arange(0, 1010, 10))/beta *
             ((np.arange(0, 1010, 10))-co2gtc[1])+tmean[1]-tmean[0],
             linewidth=2)

    ax1.axhline(tmean[1]-tmean[0], linestyle='--', color='k', linewidth=2)

    ax1.set_ylim([0, 1.5])

    ax1.set_ylabel('Global mean temperature anomaly relative to 2006-2016 [K]',
                   fontweight='bold', fontsize=20)
    ax1.set_xlabel('Cumulative total anthropogenic CO2 emissions from 2011 [GtC]', fontweight='bold', fontsize=20)
    ax1.set_xticks(np.arange(0, 1100, 100))
    ax1.set_xlim([0, 1000])
    ax1.axvline((tmean[1]-tmean[0])/g2, linestyle='--', color='r', linewidth=2)
    ax1.axvline((((alpha*dfdc(600)/beta)*co2gtc[1])+tmean[1]-tmean[0]) /
                (g2+alpha*dfdc(600)/beta),
                linestyle='-.', color='r', linewidth=2)
    ax2 = ax1.twiny()
    ax2.spines['left'].set_position(('axes', 0.0))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 1.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')
    ax2.set_xlabel('Atmospheric CO2 concentration [ppm]', fontweight='bold',
                   fontsize=20)
    ax2.set_ylim([0, 1.5])
    new_tick_locations = ppm2gtc(np.arange(400, 700, 25))

    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)

    ax2.set_xticklabels(gtc2ppm(new_tick_locations),
                        fontweight='bold', fontsize=20)

    ax1.set_yticklabels(ax1.get_yticks(), fontweight='bold', fontsize=20)
    ax1.set_xticklabels(ax1.get_xticks(), fontweight='bold', fontsize=20)
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.15, left=.05,
                        right=.97)


def boxplot(globprcnts, prcnts):
    from matplotlib import pyplot as plt
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = np.array([-2, -1, .5, 1.5, 3, 4, 5.5, 6.5, 8, 9, 11.5, 12.5, 15, 16,
                  17.5, 18.5, 21, 22, 23.5, 24.5, 27, 28, 29.5, 30.5])
    globprcnts2 = np.copy(globprcnts)
    globprcnts2[1:, 0] -= globprcnts[0, 0]
    globprcnts2 = globprcnts2[2:]
    prcnts2 = np.copy(prcnts)
    prcnts2[1:, :, 0] -= prcnts[0, :, 0]
    prcnts2 = prcnts2[2:, :, :]
    labels = ['Tglob_l', 'Tglob_u', 'Tmean_l', 'Tmean_u', 'TX_l', 'TX_u',
              'TXx_l', 'TXx_u', 'WBGT_l', 'WBGT_u', 'TX90p_l', 'TX90p_u',
              'SW_l', 'SW_u', 'LW_l', 'LW_u', 'R_l', 'R_u', 'R_trop_l',
              'R_trop_u', 'R90p_l', 'R90p_u', 'R90p_trop_l', 'R90p_trop_u']

    y = np.concatenate((globprcnts2[:, 0].flatten(order='F'),
                        prcnts2[:, :, 0].flatten(order='F')))
    yerrs = np.array([np.concatenate((globprcnts2[:, 1].flatten(order='F'),
                                      prcnts2[:, :, 1].flatten(order='F'))),
                      np.concatenate((globprcnts2[:, 2].flatten(order='F'),
                                      prcnts2[:, :, 2].flatten(order='F')))])
    ax1.errorbar(x[:10], y[:10], yerr=yerrs[:, :10], fmt='x', elinewidth=2,
                 capthick=2, capsize=5, ms=10, mew=2)

    ax2 = ax1.twinx()
    ax2.errorbar(x[10:12], y[10:12], yerr=yerrs[:, 10:12], fmt='x',
                 elinewidth=2, capthick=2, capsize=5, ms=10, mew=2)
    ax3 = ax1.twinx()
    ax3.errorbar(x[12:16], y[12:16], yerr=yerrs[:, 12:16], fmt='x',
                 elinewidth=2, capthick=2, capsize=5, ms=10, mew=2)
    ax4 = ax1.twinx()
    ax4.errorbar(x[16:20], y[16:20], yerr=yerrs[:, 16:20], fmt='x',
                 elinewidth=2, capthick=2, capsize=5, ms=10, mew=2)
    ax5 = ax1.twinx()
    ax5.errorbar(x[20:], y[20:], yerr=yerrs[:, 20:], fmt='x', elinewidth=2,
                 capthick=2, capsize=5, ms=10, mew=2)

    ax2.spines['left'].set_position(('data', 11))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 0.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax3.spines['left'].set_position(('data', 14.5))
    ax3.spines['right'].set_color('none')
    ax3.spines['bottom'].set_position(('axes', 0.0))
    ax3.spines['top'].set_color('none')
    ax3.spines['bottom'].set_smart_bounds(True)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')

    ax4.spines['left'].set_position(('data', 20.5))
    ax4.spines['right'].set_color('none')
    ax4.spines['bottom'].set_position(('axes', 0.0))
    ax4.spines['top'].set_color('none')
    ax4.spines['bottom'].set_smart_bounds(True)
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('left')
    ax4.yaxis.set_label_position('left')

    ax5.spines['left'].set_position(('data', 26.5))
    ax5.spines['right'].set_color('none')
    ax5.spines['bottom'].set_position(('axes', 0.0))
    ax5.spines['top'].set_color('none')
    ax5.spines['bottom'].set_smart_bounds(True)
    ax5.xaxis.set_ticks_position('bottom')
    ax5.yaxis.set_ticks_position('left')
    ax5.yaxis.set_label_position('left')
    ax1.set_xlim([-2.5, 31])
    ax1.set_ylabel('Temperature anomaly [K]', fontweight='bold', fontsize=20)
    ax2.set_ylabel('Percentage anomaly [%]', fontweight='bold', fontsize=20)
    ax3.set_ylabel('Flux anomaly [W/m^2]', fontweight='bold', fontsize=20)
    ax4.set_ylabel('Precipitation anomaly [mm/day]', fontweight='bold',
                   fontsize=20)
    ax5.set_ylabel('Percentage anomaly [%]', fontweight='bold', fontsize=20)
    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontweight='bold', fontsize=20)
    ax1.set_xticklabels(labels, fontweight='bold', fontsize=20,
                        rotation='vertical')
    ax2.set_yticklabels(ax2.get_yticks(), fontweight='bold', fontsize=20)
    ax3.set_yticklabels(ax3.get_yticks(), fontweight='bold', fontsize=20)
    ax4.set_yticklabels(ax4.get_yticks(), fontweight='bold', fontsize=20)
    ax5.set_yticklabels(ax5.get_yticks(), fontweight='bold', fontsize=20)
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.15, left=.05,
                        right=.99)
