# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 14:04:49 2016

@author: bakerh

NAME
    iGCM plotting functions
PURPOSE
    Defines all plotting functions for data from iGCM
"""

import numpy as np
import scipy as sp


def montecarlo(jetind, n):
    control = jetind[0, :, :]
    sig = np.zeros((2, 2, 306))
    for i in range(306):
        merge = np.concatenate((control, jetind[i+1, :]), axis=0)
        p = np.zeros((2, 2, n+1))
        for j in range(n):
            merge = np.random.permutation(merge)
            merge1 = merge[:2160, :]
            merge2 = merge[2160:, :]
            for a in range(2):
                for b in range(2):
                    p[a, b, j+1] = ((np.percentile(merge1[:, a, b], [75]) -
                                    np.percentile(merge1[:, a, b], [25])) -
                                    (np.percentile(merge2[:, a, b], [75]) -
                                    np.percentile(merge2[:, a, b], [25])))
        for a in range(2):
            for b in range(2):
                p[a, b, 0] = ((np.percentile(jetind[i+1, :, a, b], [75]) -
                              np.percentile(jetind[i+1, :, a, b], [25])) -
                              (np.percentile(control[:, a, b], [75]) -
                              np.percentile(control[:, a, b], [25])))
                if a == 0 & b == 1:
                    p[a, b, :] = p[a, b, :] * -1
                sig[a, b, i] = sp.stats.percentileofscore(p[a, b, 1:],
                                                          p[a, b, 0],
                                                          kind='mean')
    sig = (sig-50)*2/100
    sig = 1 - abs(sig)
    sig = np.reshape(sig, (2, 2, 9, 34))
    return sig


def I2Kfull(data):
    '''
    Inserts 2Kruna levels in to 2Krun levels
    '''
    data1 = np.copy(data)
    data1[35:69, :] = data[171:205, :]
    data1[69:103, :] = data[35:69, :]
    data1[103:137, :] = data[205:239, :]
    data1[137:171, :] = data[69:103, :]
    data1[171:205, :] = data[239:273, :]
    data1[205:239, :] = data[103:137, :]
    data1[239:273, :] = data[273:307, :]
    data1[273:307, :] = data[137:171, :]
    return data1


def animate(data, sigma, lat):
    """
    Plots animation of input grid of quantity at sigma and lat coords

    Parameters
    ----------
    data: array
        data being plotted
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    title: str
        optional title
    """
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    my_cmap = plt.cm.RdYlBu_r(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscale(data)
    fig = plt.figure()

    def updatefig(i):
        fig.clear()
        plot = plt.contourf(meshlat, meshsigma, data[i, :, :], ctrs,
                            cmap=newcmap, vmin=caxismin, vmax=caxismax)
        plt.gca().invert_yaxis()
        plt.yscale('linear')
        # plt.xlim([-90, 90])
        # plt.xlim([-180, 180])
        # plt.xticks(np.arange(-90, 105, 15))
        plt.title('Day ' + str((i+1)/4), y=1.08, fontsize=30)
        plt.colorbar(plot, orientation='horizontal',
                     shrink=0.5, spacing='proportional')
        plt.draw()

    anim = animation.FuncAnimation(fig, updatefig, np.ma.size(data, axis=0))
    # anim.save('test.mp4',bitrate=10000)
    return anim


def animate2(data, sigma, lat):
    """
    Plots animation of input grid of quantity at sigma and lat coords. Replots
    color bars for each plot

    Parameters
    ----------
    data: array
        data being plotted
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    title: str
        optional title
    """
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    my_cmap = plt.cm.RdYlBu_r(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    fig = plt.figure()

    def updatefig(i):
        fig.clear()
        caxismin, caxismax, ctrs = colourscale(data[i, :, :])
        plot = plt.contourf(meshlat, meshsigma, data[i, :, :], ctrs,
                            cmap=newcmap, vmin=caxismin, vmax=caxismax)
        plt.gca().invert_yaxis()
        plt.yscale('linear')
        # plt.xlim([-90, 90])
        # plt.xlim([-180, 180])
        # plt.xticks(np.arange(-90, 105, 15))
        plt.title('Run ' + str(i), y=1.08, fontsize=30)
        plt.colorbar(plot, orientation='horizontal',
                     shrink=0.5, spacing='proportional')
        plt.draw()

    anim = animation.FuncAnimation(fig, updatefig, np.ma.size(data, axis=0))
    # anim.save('test.mp4',bitrate=10000)
    return anim


def anim2(plotdata):

    from mpl_toolkits.basemap import Basemap, shiftgrid
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    lat = np.arange(90, -90-180/(np.ma.size(plotdata, axis=1)-1),
                    -180/(np.ma.size(plotdata, axis=1)-1))
    lon = np.arange(0, 360+360/(np.ma.size(plotdata, axis=2)),
                    360/(np.ma.size(plotdata, axis=2)))
    temp = np.expand_dims(plotdata[:, :, 0], axis=2)
    plotdata = np.concatenate((plotdata, temp), axis=2)
    plotdata, lon = shiftgrid(180., plotdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)
    my_cmap = plt.cm.RdYlBu_r(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscale(plotdata)
    fig = plt.figure()

    def updatefig(i):
        fig.clear()
        m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=-180, urcrnrlon=180, resolution='c')
        m.drawcoastlines()
        m.drawmapboundary()
        x, y = m(meshlon, meshlat)
        plot = m.contourf(x, y, plotdata[i, :, :], ctrs,
                          cmap=newcmap, vmin=caxismin, vmax=caxismax)
        plt.colorbar(plot, orientation='horizontal',
                     shrink=0.5, spacing='proportional')
        parallels = m.drawparallels(np.arange(-90., 91., 15.))
        meridians = m.drawmeridians(np.arange(-180., 181., 30.))
        m.drawparallels(parallels, labels=[True, True, True, True])
        m.drawmeridians(meridians, labels=[True, True, True, True])
        plt.title('Run ' + str(i), y=1.08, fontsize=30)
        plt.draw()

    anim = animation.FuncAnimation(fig, updatefig,
                                   np.ma.size(plotdata, axis=0))
    # anim.save('test.mp4',bitrate=10000)
    return anim


def barocliniczw(vT, lat):
    '''
    computes the baroclinic zone width
    from the 90th - 10th percentile latitudes
    of the v'T' dist

    Parameters:
    -----------
    vT: array
        vT
    lat: array
        latitudes

    Returns
    -------
    deltabzw: array
        difference in zone width from control
    '''
    vT850 = vT[:, 32, :]
    bzw = np.zeros((2, 307))
    for i in range(307):
        wm = np.min(vT850[i, :])
        sm = np.max(vT850[i, :])
        wam = np.argmin(vT850[i, :])
        sam = np.argmax(vT850[i, :])
        w10 = (np.abs(vT850[i, :wam]-.1*wm)).argmin()
        w90 = (np.abs(vT850[i, wam:32]-.1*wm)).argmin() + wam
        s10 = (np.abs(vT850[i, 32:sam]-.1*sm)).argmin() + 32
        s90 = (np.abs(vT850[i, sam:]-.1*sm)).argmin() + sam
        bzw[0, i] = np.abs(lat[w90] - lat[w10])
        bzw[1, i] = lat[s90] - lat[s10]
    bzw_change = np.zeros((2, 306))
    for i in range(2):
        bzw_change[i, :] = bzw[i, 1:] - bzw[i, 0]
    bzw_change = np.reshape(bzw_change, (2, 9, 34))
    return bzw_change


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
        ctrs2 = np.arange(0.1*M, 1.09*M, .1*M)
        ctrs = np.concatenate((ctrs1, ctrs2))
        caxismin = -M
        caxismax = M
    else:
        m = -m
        ctrs1 = np.arange(-m, 0, .1*m)
        ctrs2 = np.arange(0.1*m, 1.09*m, .1*m)
        ctrs = np.concatenate((ctrs1, ctrs2))
        caxismin = -m
        caxismax = m
    # function will not work if there exist no positive max or negative min
    return caxismin, caxismax, ctrs


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


def colourscaleint1(m):
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
    ctrs1 = np.arange(-m, 0, 0.1*m)
    ctrs2 = np.arange(0.1*m, 1.09*m, 0.1*m)
    ctrs = np.concatenate((ctrs1, ctrs2))
    caxismin = -m
    caxismax = m
    return caxismin, caxismax, ctrs


def corr(a, b, sigma, lat):
    '''
    Computes the weighted correlation between
    fields a and b

    Parameters
    ----------
    a: array
        field 1
    b: array
        field 2
    sigma: array
        sigma levels of fields
    lat: array
        lats of field

    Returns
    -------
    Corr: float
        correlation between a and b
    '''
    meshlat = np.zeros([np.ma.size(sigma), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.sqrt(np.cos(meshlat * np.pi/180))
    aw = a * meshlatweight
    bw = b * meshlatweight
    corr = np.sum(aw*bw) / np.sqrt(np.abs(np.sum(aw*aw)*np.sum(bw*bw)))
    return corr


def dailydist(control, jetind, gridlat):
    '''
    Plots the daily distribution of jet indices as a function
    of heating lat and level

    Parameters
    ----------
    control: array
        control run data
    jetind: array
        array of jet indices for runs
    gridlat: array
        latitudes of experiment runs
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    my_cmap = plt.cm.RdPu(np.arange(256))
    my_cmap[0:3, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    speeds = np.linspace(4.5, 23.5, 20)
    lat = np.linspace(-90.5, 90.5, 182)
    means = np.zeros((2, 2, 5, 34))
    gridlat1 = np.zeros((35))
    gridlatcopy = np.append(gridlat, gridlat[-1] + gridlat[-1] - gridlat[-2])
    gridlatcopy = np.append(gridlatcopy, gridlatcopy[-1] +
                            gridlatcopy[-1] - gridlatcopy[-2])
    for i in range(35):
        gridlat1[i] = (gridlatcopy[i] - gridlatcopy[i+1]) / 2 + gridlatcopy[i]
    gridlat1[18] = gridlat1[16] * -1
    gridlat1[15] = gridlat1[19] * -1
    diff = gridlat1[-1] - gridlat1[-2]
    gridlat1 = np.hstack((gridlat1, gridlat1[-1]+diff, gridlat1[-1]+2*diff))
    gridmean = gridlat1[-1] - (gridlat1[-1] - gridlat1[-2]) / 2
    controlmean = np.zeros((2, 2))
    for b in range(2):
        if b == 0:
            bns, lngth = speeds, np.ma.size(speeds)
            binspeeds = np.zeros((2, 5, 34, lngth))
            controlbinspeeds = np.zeros((3, lngth))
        else:
            bns, lngth = lat, np.ma.size(lat)
            binlats = np.zeros((2, 5, 34, lngth))
            controlbinlats = np.zeros((3, lngth))
        for a in range(2):
            for i in range(5):
                count = np.digitize(jetind[34*(i):34*(i+1), :, a, b], bns)
                controlcount = np.digitize(control[:, a, b], bns)
                means[a, b, i, :] = np.mean(jetind[34*(i):34*(i+1), :, a, b],
                                            axis=1)
                controlmean[a, b] = np.mean(control[:, a, b])
                for j in range(34):
                    if b == 0:
                        binspeeds[a, i, j, :] =\
                            np.bincount(count[j, :], minlength=lngth) / 21.60
                        controlbinspeeds[a+1, :] =\
                            np.bincount(controlcount, minlength=lngth) / 21.60
                    else:
                        binlats[a, i, j, :] =\
                            np.bincount(count[j, :], minlength=lngth) / 21.60
                        controlbinlats[a+1, :] =\
                            np.bincount(controlcount, minlength=lngth) / 21.60
    fig, axs = plt.subplots(4, 5, sharex='col', sharey='row', facecolor='w',
                            edgecolor='k')
    for a in range(2):
        for b in range(2):
            for c in range(5):
                if b == 0:
                    plot1 = axs[a+2*(1-b), c].pcolormesh(gridlat1, speeds, np.transpose(np.vstack((binspeeds[a, c, :, 1:],controlbinspeeds[0,1:],controlbinspeeds[a+1,1:]))), cmap=newcmap, vmax=50)
                    axs[a+2*(1-b), c].scatter(gridlat, means[a, b, c, :], marker='x', color='Blue', linewidth=1, s=20)
                    axs[a+2*(1-b), c].scatter(gridmean, controlmean[a, b], marker='x', color='Blue', linewidth=1, s=20)
                    if a == 1:
                        axs[a+2*(1-b), c].set_xlabel('Heating latitude (deg)')
                    if c == 0:
                        axs[a+2*(1-b), c].set_ylabel('Jet speed (m/s)')
                else:
                    plot2 = axs[a+2*(1-b), c].pcolormesh(gridlat1, lat, np.transpose(np.vstack((binlats[a, c, :, 1:],controlbinlats[0,1:],controlbinlats[a+1,1:]))), cmap=newcmap, vmax=25)
                    axs[a+2*(1-b), c].scatter(gridlat, means[a, b, c, :], marker='x', color='Blue', linewidth=1, s=20)
                    axs[a+2*(1-b), c].scatter(gridmean, controlmean[a, b], marker='x', color='Blue', linewidth=1, s=20)                    
                    if c == 0:
                        axs[a+2*(1-b), c].set_ylabel('Jet latitude (deg)')
                axs[1, 0].set_ylim(40, 80)
                axs[0, 0].set_ylim(-70, -30)
                axs[2, 0].set_ylim(speeds[0], speeds[-1])
                axs[3, 0].set_ylim(speeds[0], speeds[-1])
                axs[0, c].set_xlim(gridlat1[0], gridlat1[-1])
                axs[0, c].set_xticks(np.arange(-60, 90, 30))
                plt.subplots_adjust(hspace=.11, wspace=.1)
    cbar_ax = fig.add_axes([0.925, 0.1, 0.01, 0.38])
    fig.colorbar(plot1, cax=cbar_ax, spacing='proportional',
                 label="relative frequency (%)", extend='max')
    cbar_ax1 = fig.add_axes([0.925, 0.51, 0.01, 0.38])
    fig.colorbar(plot2, cax=cbar_ax1, spacing='proportional',
                 label="relative frequency (%)", extend='max')
    cols = ['1000hPa', '800hPa', '600hPa', '400hPa', '200hPa']
    rows = ['Winter jet', 'Summer jet', 'Winter jet', 'Summer jet']

    pad = 5  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    size='large', ha='center', va='baseline')
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation=90)


def dailydistp(control, jetind, gridlat):
    '''
    Plots the daily distribution of jet indices as a function
    of heating lat and level in percentiles

    Parameters
    ----------
    control: array
        control run data
    jetindices: array
        array of jet indices for runs
    gridlat: array
        latitudes of experiment runs
    '''
    import matplotlib.pyplot as plt
    controlp = np.zeros((2, 2, 5))
    p = np.zeros((2, 2, 9, 34, 5))
    ptls = [10, 30, 50, 70, 90]
    for a in range(2):
        for b in range(2):
            controlp[a, b, :] = np.percentile(control[:, a, b], ptls)
            for c in range(9):
                p[a, b, c, :, :] = np.transpose(np.percentile
                                                (jetind[34*(c):34*(c+1),
                                                        :, a, b], ptls,
                                                 axis=1))

    fig, axs = plt.subplots(4, 9, sharex='col', sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    for a in range(2):
        for b in range(2):
            for c in range(9):
                axs[2*a+1-b, c].plot(gridlat, p[a, b, c, :, 0], color='r',
                                     linewidth=2)
                axs[2*a+1-b, c].plot(gridlat, p[a, b, c, :, 1], color='r',
                                     linewidth=2, linestyle='dashed')
                axs[2*a+1-b, c].plot(gridlat, p[a, b, c, :, 3], color='r',
                                     linewidth=2, linestyle='dashed')
                axs[2*a+1-b, c].plot(gridlat, p[a, b, c, :, 4], color='r',
                                     linewidth=2)
                axs[2*a+1-b, c].plot(gridlat, p[a, b, c, :, 2], color='r',
                                     linewidth=2)
                axs[2*a+1-b, c].axhline(controlp[a, b, 2], color='black',
                                        linewidth=2)
                axs[2*a+1-b, c].axhline(controlp[a, b, 1], color='black',
                                        linewidth=2, linestyle='dashed')
                axs[2*a+1-b, c].axhline(controlp[a, b, 4], color='black',
                                        linewidth=2)
                axs[2*a+1-b, c].axhline(controlp[a, b, 3], color='black',
                                        linewidth=2, linestyle='dashed')
                axs[2*a+1-b, c].axhline(controlp[a, b, 0], color='black',
                                        linewidth=2)
                axs[2*a+1-b, c].axvline(controlp[a, 1, 0], color='blue',
                                        linewidth=2)
                if b == 0:
                    if c == 0:
                        axs[2*a+1-b, c].set_ylabel('Jet speed (m/s)',
                                                   fontweight='bold')
                else:
                    if c == 0:
                        axs[2*a+1-b, c].set_ylabel('Jet latitude (deg)',
                                                   fontweight='bold')
                axs[3, 4].set_xlabel('Heating latitude (deg)',
                                     fontweight='bold', fontsize=16)
                axs[2, 0].set_ylim(43, 79)
                axs[0, 0].set_ylim(-68, -32)
                axs[1, 0].set_ylim(9, 19)
                axs[3, 0].set_ylim(5, 15)
                axs[0, c].set_xlim(gridlat[0], gridlat[-1])
                axs[0, c].set_xticks(np.arange(-60, 90, 30))
                plt.subplots_adjust(hspace=.11, wspace=.1)
                for axis in ['top', 'bottom', 'left', 'right']:
                    axs[2*a+1-b, c].spines[axis].set_linewidth(2)
                axs[2*a+1-b, c].xaxis.set_tick_params(width=2)
                axs[2*a+1-b, c].yaxis.set_tick_params(width=2)
    cols = ['1000hPa', '900hPa', '800hPa', '700hPa', '600hPa', '500hPa',
            '400hPa', '300hPa', '200hPa']
    rows = ['Winter jet', '', 'Summer jet']

    pad = 5  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    size='large', ha='center', va='baseline',
                    fontweight='bold')
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, -100),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center', rotation=90,
                    fontweight='bold', fontsize=16)


def div(field, sigma, lat):
    """
    Takes the meridional divergence of input data

    Parameters
    ----------
    field: array
        data field
    sigma: array
        sigma levels of field
    lat: array
        lat of field

    Returns
    -------
    divF: array
        divergence of field
    """
    import numpy as np
    #  dFsig = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    dFlat = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    dlat = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    #  dsig = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    '''
    for i in range(np.ma.size(sigma)-1):
        dFsig[i, :] = field[i+1, :] - field[i, :]
    dFsig[-1, :] = dFsig[-2, :]
    for i in range(np.ma.size(sigma)-1):
        dsig[i, :] = meshsigma[i+1, :] - meshsigma[i, :]
    dsig[-1, :] = dsig[-2, :]
    '''
    fieldw = field * np.cos(meshlat*np.pi/180)
    for i in range(np.ma.size(lat)-1):
        dFlat[:, i] = fieldw[:, i+1] - fieldw[:, i]
    dFlat[:, -1] = dFlat[:, -2]
    for i in range(np.ma.size(lat)-1):
        dlat[:, i] = meshlat[:, i+1] - meshlat[:, i]
    dlat[:, -1] = dlat[:, -2]
    divF = dFlat * 180 / (dlat * np.cos(meshlat*np.pi/180) * 6371000 * np.pi)
    return divF


def efold(data):
    '''
    Computes the efolding time of input data

    Parameters
    ----------
    data: array
        data to compute efolding time for

    Returns
    -------
    efold: array
        efolding times at each spatial location
    '''
    import scipy as sp
    from scipy import optimize
    field = data.copy()
    n = 60
    auto = np.zeros((n))
    for i in range(n):
        auto[i] = np.corrcoef(field[n-1:-1], field[n-1-i:-i-1])[0, 1]

    def func(x, a):
        return np.exp(-x/a)

    tau = sp.optimize.curve_fit(func, np.linspace(0, n-1, n), auto)[0]
    return auto, tau


def eof_analysis(data, sigma, lat, u):
    '''
    Amalgamates all EOF functions into one process

    Parameters
    ----------
    data: array
        data to determine EOFs for. time x nsigma x nlat
    sigma: array
        sigma levels of data
    lat: array
        latitude of data
    u: array
        results from experiment to project on to EOF

    Returns
    -------
    eofs: array
        First 3 EOFs. 3 x nsigma x nlat
    var: array
        amount of variance explained by each EOF. 20 x 2 array:
        var(:, 1) = the 1st 20 eigenvalues of the covariance matrix
         ar(:, 2) = the 1st 20 explained variances (in %).
    response: array
        response in mean zonal wind projected onto the first 3 EOFs
    '''
    nth = 3
    field = np.copy(data)
    sig = np.copy(sigma)
    lati = np.copy(lat)
    eofs, pcs, var = eof_svd(field, sig, lati)
    eofreg = np.zeros((nth, 37, np.ma.size(lati)))
    for i in range(nth):
        eofreg[i, :, :] = eof_regress(pcs, i+1, field)
    response = np.zeros((nth, 9, 34))
    eof = np.zeros((nth, 37, np.ma.size(lati)))
    for i in range(nth):
        response[i, :], eof[i, :] = eof_response(u, eofreg[i, :], sig, lati)
    return response, eof, var


def eof_svd(data, sigma, lat):
    '''
    EOFs of data

    Parameters
    ----------
    data: array
        data to determine EOFs for. time x nsigma x nlat
    sigma: array
        sigma levels of data
    lat: array
        latitude of data

    Returns
    -------
    eofs: array
        First 20 EOFs. 20 x nsigma x nlat
    var: array
        amount of variance explained by each EOF. 5 x 2 array:
        var(:, 1) = the 1st 5 eigenvalues of the covariance matrix
         ar(:, 2) = the 1st 5 explained variances (in %).
    pc: array
        principal components. time x 5
    '''
    field = data.copy()
    from numpy import linalg
    neof = 20
    dim = field.shape
    field = eofweight(field, sigma, lat)
    aver = np.mean(field, axis=0)
    for i in range(dim[0]):
        field[i, :] = field[i, :] - aver
    field = np.reshape(field, (dim[0], dim[1] * dim[2]), order='F')

    if dim[0] > dim[1] * dim[2]:
        field = np.transpose(field)

    u, s, v = linalg.svd(field)
    v = np.transpose(v)
    s = s * s
    eigs = np.copy(s)
    s = s / np.sum(s)
    if dim[0] < dim[1] * dim[2]:
        u, v = v, u

    eofs = u[:, :neof]
    pcs = v[:, :neof]
    var = np.zeros([neof, 2])
    var[:, 0] = eigs[:neof]
    var[:, 1] = s[:neof] * 100
    eofs = np.transpose(eofs)
    eofs = np.reshape(eofs, (neof, dim[1], dim[2]), order='F')
    return eofs, pcs, var


def eof_regress(pcs, neof, field):
    '''
    Regresses original field on nth pcs time series

    Parameters
    ----------
    field: array
        data to regress. time x nsigma x nlat
    neof: int
        eof number to calculate
    pcs: array
        first principal component time series

    Returns
    -------
    eofn: array
        nth eof from field regressed onto pcsn
    '''
    from scipy import stats
    pcsn = pcs[:, neof-1]
    pcsn = (pcsn - np.mean(pcsn)) / np.std(pcsn)
    eofn = np.zeros((np.ma.size(field, axis=1),
                    np.ma.size(field, axis=2)))
    for a in range(np.ma.size(field, axis=1)):
        for b in range(np.ma.size(field, axis=2)):
            eofn[a, b] = stats.linregress(pcsn, field[:, a, b])[0]
    return eofn  # this is not normalized!


def eof_response(u, eof1, sigma, lat):
    '''
    Projects response of zonal wind onto eof1

    Parameters
    ----------
    u: array
        data to project on to eof1
    eof1: array
        first eof
    sigma: array
        sigma levels of data
    lat: array
        latitude of data

    Returns
    -------
    response: array
        response projected on eof1
    gridlat: array
        lat points used in experiment
    gridsigma: array
        sigma points used in experiment
    '''

    dim = eof1.shape
    du = np.zeros((306, 37, dim[1]))
    for i in range(306):
        du[i, :] = u[i+1, :, :] - u[0, :, :]
    du = eofweight(du, sigma, lat)
    eof1w = eofweight(eof1, sigma, lat)  # units here?? what if sig 1000 times
    eof = eof1 / ((np.sum(eof1w * eof1w)))  # normalized eof TAKEN SQRTS OUT _
    eof1w = eof1w / np.sqrt((np.sum(eof1w * eof1w)))  # norm and weighted for use
    response = np.zeros((306))
    for i in range(306):
        response[i] = np.sum(du[i, :, :] * eof1w)
    response = np.reshape(response, (9, 34))
    return response, eof


def eofweight(unweighted, sigma, lat):
    '''
    Outputs weighted data

    Parameters
    ----------
    unweighted: array
        unweighted data
    lat: array
        latitudes
    sigma: array
        sigma values

    Outputs
    -------
    weighted: array
        weighted data
    '''
    sig = sigma.copy()
    meshlat = np.zeros([np.ma.size(sigma), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.sqrt(np.cos(meshlat * np.pi/180))
    sig = np.insert(sig, 0, 0)
    sig = np.append(sig, 1)
    deltasigma = np.zeros(np.ma.size(sigma))
    for i in range(np.ma.size(sigma)):
        deltasigma[i] = (sig[i+2] - sig[i]) / 2
    meshsigma = np.zeros([np.ma.size(lat), np.ma.size(sigma)])
    meshsigma[:, :] = deltasigma
    meshsigma = np.sqrt(np.transpose(meshsigma))
    weighted = unweighted * meshlatweight * meshsigma
    return weighted


def grids(sigma, lat):
    '''
    Outputs heating grid for experiment

    Parameters
    ----------
    lat: array
        latitudes
    sigma: array
        sigma values

    Outputs
    -------
    gridlat: array
        lat array
    gridsigma: array
        sigma array
    '''
    gridlat = np.zeros((14))
    gridlat[0], gridlat[1], gridlat[2] = lat[1], lat[6], lat[11]
    gridlat[3], gridlat[4], gridlat[5] = lat[16], lat[21], lat[26]
    gridlat[6], gridlat[7], gridlat[8] = lat[31], lat[32], lat[37]
    gridlat[9], gridlat[10], gridlat[11] = lat[42], lat[47], lat[52]
    gridlat[12], gridlat[13] = lat[57], lat[62]
    gridsigma = np.zeros((5))
    gridsigma[0], gridsigma[1], gridsigma[2] = sigma[12], sigma[18], sigma[24]
    gridsigma[3], gridsigma[4] = sigma[30], sigma[36]
    return gridsigma, gridlat


def heatingarray(paras, sigma, lat, q0, sig_x, sig_y):
    """
    Outputs array contatining heating blobs

    Parameters
    ----------
    paras: array
        array containing positions of blobs
    lat: array
        latitude values for data
    sigma: array
        sigma values for data

    Outputs
    -------
    heating: array
        array of heating profiles for each run
    """
    import numpy as np
    sig_xdeg = sig_x * 180/np.pi
    heating = np.zeros((np.ma.size(paras, axis=0)-1, 30, 64))
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    for i in range(np.ma.size(paras, axis=0)-1):
        heating[i, :, :] = (q0 *
                            np.exp(-(((meshlat - 180*float(paras[i+1, 1]) /
                                   np.pi)**2) / (2*(sig_xdeg**2)) +
                                   ((meshsigma - float(paras[i+1, 2]))**2) /
                                   (2*(sig_y**2))))/86400)
    return heating


def jetindicessummer(u, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    iGCM output

    Parameters
    ----------
    u: array
        data to determine jet indices for
    lat: array
        latitude values for data

    Returns
    -------
    jetspeed: float
        max jet speed
    jetlatitude: float
        position of maximum
    """
    umean = u[:, :, 32:]
    speeds = np.zeros((np.ma.size(umean, axis=0), np.ma.size(umean, axis=1)))
    lats = np.zeros((np.ma.size(umean, axis=0), np.ma.size(umean, axis=1)))
    for l in range(np.ma.size(umean, axis=1)):
        for t in range(np.ma.size(umean, axis=0)):
            speeds[t, l] = np.amax(umean[t, l, :])
            lats[t, l] = lat[np.argmax(umean[t, l, 7:])+39]
    return speeds, lats


def jetindiceswinter(u, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    iGCM output

    Parameters
    ----------
    u: array
        data to determine jet indices for
    lat: array
        latitude values for data

    Returns
    -------
    jetspeed: float
        max jet speed
    jetlatitude: float
        position of maximum
    """
    umean = u[:, :, 0:32]
    speeds = np.zeros((np.ma.size(umean, axis=0), np.ma.size(umean, axis=1)))
    lats = np.zeros((np.ma.size(umean, axis=0), np.ma.size(umean, axis=1)))
    for l in range(np.ma.size(umean, axis=1)):
        for t in range(np.ma.size(umean, axis=0)):
            speeds[t, l] = np.amax(umean[t, l, :])
            lats[t, l] = lat[np.argmax(umean[t, l, :])]
    return speeds, lats


def lapserate(t, z, sigma, lat):
    """
    Produces plot of lapse rate of T data

    Parameters
    ----------
    t: array
        temperature data field
    z: array
        geopotential height of field
    sigma: array
        sigma levels of field
    lat: array
        lat of field
    """
    import numpy as np
    dT = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    dz = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    for i in range(np.ma.size(sigma, axis=0)-1):
        dT[i, :] = t[i+1, :] - t[i, :]
    for i in range(np.ma.size(sigma, axis=0)-1):
        dz[i, :] = z[i+1, :] - z[i, :]
    lapse = -1000 * dT[0:-1] / dz[0:-1]
    # zonalplot(lapse, sigma[0:-1], lat, 'Lapse rate')
    return lapse


def paper1_q0plot(d_q0):
    import matplotlib.pyplot as plt
    xax = np.array([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0])
    fig, axs = plt.subplots(1, 2, facecolor='w', edgecolor='k', linewidth=2)
    plt.suptitle('Linearity of response with respect to magnitude of heating',
                 y=.95, fontsize=30, fontweight='bold')
    axs[0].plot(xax, d_q0[0, :, 0, 0], ms=8, color='r', marker='o',
                linewidth=2, label='Winter 1')

    axs[0].plot(xax, d_q0[1, :, 0, 0], ms=8, color='b', marker='o',
                linewidth=2, label='Winter 2')

    axs[0].plot(xax, d_q0[2, :, 0, 0], ms=8, color='g', marker='o',
                linewidth=2, label='Winter 3')

    axs[0].plot(xax, d_q0[3, :, 0, 0], ms=8, color='k', marker='o',
                linewidth=2, label='Winter 4')

    axs[0].plot(xax, d_q0[0, :, 1, 0], ms=8, mew=2, ls='--', color='r',
                marker='x', linewidth=2, label='Summer 1')

    axs[0].plot(xax, d_q0[1, :, 1, 0], ms=8, mew=2, ls='--', color='b',
                marker='x', linewidth=2, label='Summer 2')

    axs[0].plot(xax, d_q0[2, :, 1, 0], ms=8, mew=2, ls='--', color='g',
                marker='x', linewidth=2, label='Summer 3')

    axs[0].plot(xax, d_q0[3, :, 1, 0], ms=8, mew=2, ls='--', color='k',
                marker='x', linewidth=2, label='Summer 4')
    axs[0].set_xlim(0, 10.1)
    axs[0].set_ylim(-7.5, 4)
    axs[0].set_xticks([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0])
    axs[0].set_xticklabels(axs[0].get_xticks(), fontweight='bold', fontsize=20)
    axs[0].set_yticklabels(axs[0].get_yticks(), fontweight='bold', fontsize=20)
    axs[0].set_xlabel('Heating [K/day]', weight='bold', size=20)
    axs[0].set_ylabel('Change in jet speed [m/s]', weight='bold', size=20)
    handles, labels = axs[0].get_legend_handles_labels()
    axs[0].legend(handles, labels, loc=2)

    axs[1].plot(xax, d_q0[0, :, 0, 1], ms=8, color='r', marker='o',
                linewidth=2, label='Winter 1')

    axs[1].plot(xax, d_q0[1, :, 0, 1], ms=8, color='b', marker='o',
                linewidth=2, label='Winter 2')

    axs[1].plot(xax, d_q0[2, :, 0, 1], ms=8, color='g', marker='o',
                linewidth=2, label='Winter 3')

    axs[1].plot(xax, d_q0[3, :, 0, 1], ms=8, color='k', marker='o',
                linewidth=2, label='Winter 4')

    axs[1].plot(xax, d_q0[0, :, 1, 1], ms=8, mew=2, ls='--', color='r',
                marker='x', linewidth=2, label='Summer 1')

    axs[1].plot(xax, d_q0[1, :, 1, 1], ms=8, mew=2, ls='--', color='b',
                marker='x', linewidth=2, label='Summer 2')

    axs[1].plot(xax, d_q0[2, :, 1, 1], ms=8, mew=2, ls='--', color='g',
                marker='x', linewidth=2, label='Summer 3')

    axs[1].plot(xax, d_q0[3, :, 1, 1], ms=8, mew=2, ls='--', color='k',
                marker='x', linewidth=2, label='Summer 4')
    axs[1].set_xlim(0, 10.1)
    axs[1].set_ylim(-12.5, 23)
    axs[1].set_xticks([0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 7.5, 10.0])
    axs[1].set_xticklabels(axs[1].get_xticks(), fontweight='bold', fontsize=20)
    axs[1].set_yticklabels(axs[1].get_yticks(), fontweight='bold', fontsize=20)
    axs[1].set_xlabel('Heating [K/day]', weight='bold', size=20)
    axs[1].set_ylabel('Poleward shift of jet [deg]', weight='bold', size=20)
    handles, labels = axs[1].get_legend_handles_labels()
    axs[1].legend(handles, labels, loc=2)


def paper1_qaddplot(d_qaddlin, d_qadd, jetind_solo):

    import matplotlib.pyplot as plt
    a = 7
    fig, axs = plt.subplots(1, 2, facecolor='w', edgecolor='k', linewidth=2)
    plt.suptitle('Linearity of combined responses',
                 y=.95, fontsize=30, fontweight='bold')
    axs[0].plot([-7, 17], [-7, 17], c='k', linewidth=2, alpha=0.5)    
    axs[0].scatter(d_qaddlin[:, 0, 0], d_qadd[:, 0, 0], s=15, color='b',
                   marker='o', linewidth=a, label='Winter')
    axs[0].scatter(jetind_solo[:, 0, 0], jetind_solo[:, 0, 0], s=100, facecolors='none', edgecolors='b',
                   marker='o', linewidth=1)

    axs[0].scatter(d_qaddlin[:, 1, 0], d_qadd[:, 1, 0], s=15, color='r',
                   marker='o', linewidth=a, label='Summer')
    axs[0].scatter(jetind_solo[:, 1, 0], jetind_solo[:, 1, 0], s=100, facecolors='none', edgecolors='r',
                   marker='o', linewidth=1)

    axs[0].set_xlim(-4, 1.5)
    axs[0].set_ylim(-4, 1.5)  
    axs[0].set_xticklabels(axs[0].get_xticks(), fontweight='bold', fontsize=20)
    axs[0].set_yticklabels(axs[0].get_yticks(), fontweight='bold', fontsize=20)
    axs[0].set_xlabel('Change in jet speed in sum of individual simulations [m/s]', weight='bold', size=20)
    axs[0].set_ylabel('Change in jet speed in combined forcing runs [m/s]', weight='bold', size=20)  
    handles, labels = axs[0].get_legend_handles_labels()
    axs[0].legend(handles, labels, loc=2, scatterpoints=1)
    axs[0].axhline(0, ls='--', color='gray')
    axs[0].axvline(0, ls='--', color='gray')
    labels1 = ['1', '2', '3', '4', '1 ', '2 ', '3 ', '4 ', '1+2', '1+3', '1+4',
               '2+3', '2+4', '3+4', '1+2+3', '1+2+4', '2+3+4', '1+3+4',
               '1+2+3+4', '1+2 ', '1+3 ', '1+4 ', '2+3 ', '2+4 ',
               '3+4 ', '1+2+3 ', '1+2+4 ', '2+3+4 ', '1+3+4 ', '1+2+3+4 ']

    '''
    for x, y, s in zip(np.concatenate((jetind_solo[:, 0, 0], jetind_solo[:, 1, 0], d_qaddlin[:, 0, 0],
                                       d_qaddlin[:, 1, 0])),
                       np.concatenate((jetind_solo[:, 0, 0], jetind_solo[:, 1, 0], d_qadd[:, 0, 0], d_qadd[:, 1, 0])),
                       labels1):
        axs[0].annotate(s, xy=(x, y), xytext=(x-.1*x, y+.1*y), fontweight='bold', fontsize=15,
                        arrowprops=dict(arrowstyle="->", color='k', lw=1)).draggable()
    '''
       
    axs[1].plot([-15, 17], [-15, 17], c='k', linewidth=2, alpha=0.5)
    axs[1].scatter(d_qaddlin[:, 0, 1], d_qadd[:, 0, 1], s=15, color='b',
                   marker='o', linewidth=a, label='Winter')

    axs[1].scatter(jetind_solo[:, 0, 1], jetind_solo[:, 0, 1], s=100, facecolors='none', edgecolors='b',
                   marker='o', linewidth=1)

    axs[1].scatter(d_qaddlin[:, 1, 1], d_qadd[:, 1, 1], s=15, color='r',
                   marker='o', linewidth=a, label='Summer')

    axs[1].scatter(jetind_solo[:, 1, 1], jetind_solo[:, 1, 1], s=100, facecolors='none', edgecolors='r',
                   marker='o', linewidth=1)

    axs[1].set_xlim(-5, 17)
    axs[1].set_ylim(-5, 17)   
    axs[1].set_xticklabels(axs[1].get_xticks(), fontweight='bold', fontsize=20)
    axs[1].set_yticklabels(axs[1].get_yticks(), fontweight='bold', fontsize=20)
    axs[1].set_xlabel('Poleward shift of jet in sum of individual simulations [deg]', weight='bold', size=20)
    axs[1].set_ylabel('Poleward shift in combined forcing runs [deg]', weight='bold', size=20)
    handles, labels = axs[1].get_legend_handles_labels()
    axs[1].legend(handles, labels, loc=2, scatterpoints=1)
    axs[1].axhline(0, ls='--', color='gray')
    axs[1].axvline(0, ls='--', color='gray')


def surftemp(temp, lat, title=''):
    """
    Outputs surface temperature profile

    Parameters
    ----------
    temp: array
        temp field containing all temps
    lat: array
        latitude values for data
    """
    import matplotlib.pyplot as plt
    plt.figure()
    surftemp = temp-273.15
    plt.plot(lat, surftemp)
    plt.xlim([-90, 90])
    plt.xticks(np.arange(-90, 105, 15))
    plt.title(title, y=1.08, fontsize=30)
    plt.show()


def timeplot(series, title='', colour='blue'):
    '''
    timeplot takes a times series and outputs a plot of the timeseries from
    the start date

    Parameters
    ----------
    series : array
    colour : str
        colour of data line

    Returns
    -------
    plot : fig
        plot of time series
    '''
    import matplotlib.pyplot as plt
    plt.figure()
    x = np.arange(1, np.ma.size(series)+1)
    plt.plot(x, series, color=colour)
    # plt.axhline(0, color='black')
    # plt.axvline(84, color='black')
    # plt.axvline(96, color='black')
    plt.title(title, y=1.08, fontsize=30)
    plt.show()


def zonalplot(plotdata, sigma, lat, maxx=1, title=''):
    """
    Plots input grid of quantity at sigma and lat coords

    Parameters
    ----------
    plotdata: array
        data being plotted
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    title: str
        optional title
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.figure(facecolor='white')
    #  plt.suptitle("Tropical + polar heating", size=36, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(-maxx, maxx, 17)
    plot = plt.contourf(meshlat, meshsigma, plotdata, ctrs, extend='both',
                        cmap=newcmap, vmin=-maxx, vmax=maxx)
    c = plt.colorbar(plot, orientation='horizontal',
                     shrink=1, spacing='proportional', aspect=50,
                     ticks=np.linspace(-maxx, maxx, 9))
    c.set_label('Temperature gradient difference (K/deg)', fontsize='20')
    c.ax.set_xticklabels(np.linspace(-maxx, maxx, 9), fontsize='20')
    plt.gca().invert_yaxis()
    plt.yscale('linear')
    ax = plt.gca()
    ax.set_xlim([-90, 90])
    plt.xlabel('Latitude (deg)', fontsize=20)
    plt.ylabel('Pressure (hPa)', fontsize=20)
    plt.xticks(np.arange(-90, 105, 15), fontsize='20')
    plt.yticks(fontsize='20')
    plt.title(title, y=1.05, fontsize=30)
    plt.show()


def zonalplotisen(plotdata, potemp, sigma, lat, title=''):
    """
    Plots input grid of quantity at sigma and lat coords

    Parameters
    ----------
    plotdata: array
        data being plotted
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    title: str
        optional title
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.figure()
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscaleint(10)
    plot = plt.contourf(meshlat, meshsigma, plotdata, ctrs, extend='both',
                        cmap=newcmap, vmin=caxismin, vmax=caxismax)
    c = plt.colorbar(plot, orientation='horizontal',
                     shrink=0.5, spacing='proportional',
                     ticks=np.linspace(-10, 10, 9))
    c.set_label('Temperature [C]', fontweight='bold', fontsize='24')
    c.ax.set_xticklabels(np.linspace(-10, 10, 9), fontweight='bold',
                         fontsize='20')
    ctrs2 = [270, 280, 290, 300, 310, 320, 330, 340, 350, 400, 500]
    plot1 = plt.contour(meshlat, meshsigma, potemp, ctrs2,
                        colors='k', linewidths=1.5)
    plt.clabel(plot1, inline=True, inline_spacing=0, fontsize=10, fmt='%.0f')
    plt.gca().invert_yaxis()
    plt.yscale('linear')
    plt.xlim([-90, 90])
    plt.xlabel('Latitude [deg]', fontsize=20, fontweight='bold')
    plt.ylabel('Pressure [hPa]', fontsize=20, fontweight='bold')
    plt.xticks(np.arange(-90, 105, 15), weight='bold', fontsize='20')
    plt.yticks(weight='bold', fontsize='20')
    plt.title(title, y=1.05, fontsize=30, fontweight='bold')
    plt.show()


def zonalplotmark(plotdata, p, sigma, lat, control, csigma, clat, title='',
                  clabel=''):
    """
    Plots input grid of quantity at sigma and lat coords
    Marks data points with crosses

    Parameters
    ----------
    plotdata: array
        data being plotted
    p: array
        significance of runs
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    control: array
        control data
    csigma: array
        sigma levels of control
    clat: array
        latitudes of control
    title: str
        optional title
    clabel: str
        optional colourbar caption
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    maxes = np.max(plotdata), np.min(plotdata)
    maxes1 = np.max(control), np.min(control)
    m = float(input('Enter colourbar max given %f, %f: ' % (maxes)))
    m1 = float(input('Enter contour max given %f, %f: ' % (maxes1)))
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    plt.figure()
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscaleint(m)
    plot = plt.contourf(meshlat, meshsigma, plotdata, ctrs, extend='both',
                        cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional', format='%.1f',
                 label=clabel)
    cmi, cma, cctrs = colourscaleint(m1)
    plot1 = plt.contour(cmeshlat, cmeshsigma, control, cctrs,
                        colors='k', linewidths=1.5)
    plt.clabel(plot1, inline=True, inline_spacing=-3, fontsize=10, fmt='%.0f')
    plt.yscale('linear')
    plt.gca().invert_yaxis()
    plt.xlim([-90, 90])
    plt.ylim([1, 0])
    plt.xticks(np.arange(-90, 105, 15))
    plt.title(title, y=1.08, fontsize=30)
    # crosses
    meshsigma1 = np.copy(meshsigma)
    for a in range(9):
        for b in range(34):
            if p[a, b] <= 0.005:
                meshsigma[a, b] = np.nan
    plt.scatter(meshlat, meshsigma, marker='o', color='gray', s=0.5)
    for a in range(9):
        for b in range(34):
            if p[a, b] > 0.005:
                meshsigma1[a, b] = np.nan
    plt.scatter(meshlat, meshsigma1, marker='o', color='black', s=5)
    plt.show()


def zonalusens(data, p, sigma, lat, control, csigma, clat):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    plt.suptitle("Jet sensitivity", size=36, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscaleint(5)
    caxismin1, caxismax1, ctr1 = colourscaleint(10.0)
    cmi, cma, cctrs = colourscaleint(50)

    plot = axs[0, 0].contourf(meshlat, meshsigma, data[0, 1, :], ctr1,
                              extend='both',
                              cmap=newcmap, vmin=caxismin1, vmax=caxismax1)
    plot1 = axs[0, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot = axs[0, 1].contourf(meshlat, meshsigma, data[0, 0, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plot1 = axs[0, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot = axs[1, 0].contourf(meshlat, meshsigma, data[1, 1, :], ctr1,
                              extend='both', cmap=newcmap, vmin=caxismin1,
                              vmax=caxismax1)
    plot1 = axs[1, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot2 = axs[1, 1].contourf(meshlat, meshsigma, data[1, 0, :], ctrs,
                               extend='both', cmap=newcmap, vmin=caxismin,
                               vmax=caxismax)
    plot1 = axs[1, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    for i in range(2):
        for j in range(2):
            axs[i, j].invert_yaxis()
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(),
                                      fontweight='bold', fontsize=20)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(),
                                      fontweight='bold', fontsize=20)
            meshsigma2 = np.copy(meshsigma)
            meshsigma1 = np.copy(meshsigma)
            # axs[i, j].contour(meshlat, meshsigma, p[i, 1-j, :], [0.005],
            # linewidths=2)
            for a in range(9):
                for b in range(34):
                    if p[i, 1-j, a, b] <= 0.005:
                        meshsigma2[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma2, marker='o', color='gray',
                              s=0.5)
            for a in range(9):
                for b in range(34):
                    if p[i, 1-j, a, b] > 0.005:
                        meshsigma1[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma1, marker='o', color='black',
                              s=5)
    plt.subplots_adjust(hspace=0.079, wspace=.05)
    cbar_ax = fig.add_axes([0.125, 0.05, 0.38, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin1, caxismax1, 9))
    b.set_label(label='Poleward jet latitude shift [deg]', size=20,
                weight='bold')
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cbar_ax1 = fig.add_axes([0.5225, 0.05, 0.38, 0.015])
    a = fig.colorbar(plot2, cax=cbar_ax1, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin, caxismax, 9))
    a.set_label(label='Change in amplitude of jet speed [m/s]', size=20, weight='bold')
    cl = plt.getp(cbar_ax1, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cols = ['Latitude', 'Magnitude']
    rows = ['Winter response', 'Summer response']
    pad = 10  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline',
                    fontweight='bold', fontsize=24)
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90,
                    fontweight='bold', fontsize=24)

    plt.show()


def zonalduvsens(data1, p, sigma, lat, control1, heat, csigma, clat):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    plt.suptitle("Sensitivity of u'v' divergence", size=30, fontweight='bold',
                 y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscaleint(20)
    caxismin1, caxismax1, ctr1 = colourscaleint(10)
    cctrs = [-70, -60, -50, -40, -30, -20, -10, 10, 20, 30, 40, 50]
    data = np.copy(data1)
    data[:, 0, :] = data[:, 0, :] * 1e6
    control = control1 * 1e6

    plot = axs[0, 0].contourf(meshlat, meshsigma, data[0, 1, :], ctr1,
                              extend='both',
                              cmap=newcmap, vmin=caxismin1, vmax=caxismax1)
    plot1 = axs[0, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot = axs[0, 1].contourf(meshlat, meshsigma, data[0, 0, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plot1 = axs[0, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot = axs[1, 0].contourf(meshlat, meshsigma, data[1, 1, :], ctr1,
                              extend='both', cmap=newcmap, vmin=caxismin1,
                              vmax=caxismax1)
    plot1 = axs[1, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot2 = axs[1, 1].contourf(meshlat, meshsigma, data[1, 0, :], ctrs,
                               extend='both', cmap=newcmap, vmin=caxismin,
                               vmax=caxismax)
    plot1 = axs[1, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    for i in range(2):
        for j in range(2):
            axs[i, j].invert_yaxis()
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(),
                                      fontweight='bold', fontsize=20)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(),
                                      fontweight='bold', fontsize=20)
            meshsigma2 = np.copy(meshsigma)
            meshsigma1 = np.copy(meshsigma)
            for a in range(9):
                for b in range(34):
                    if p[i, 1-j, a, b] <= 0.005:
                        meshsigma2[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma2, marker='o', color='gray',
                              s=0.5)
            for a in range(9):
                for b in range(34):
                    if p[i, 1-j, a, b] > 0.005:
                        meshsigma1[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma1, marker='o', color='black',
                              s=5)
    plt.subplots_adjust(hspace=0.079, wspace=.05)
    cbar_ax = fig.add_axes([0.125, 0.05, 0.38, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin1, caxismax1, 9))
    b.set_label(label="Poleward shift in divergence of u'v' [deg]", size=20,
                weight='bold')
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cbar_ax1 = fig.add_axes([0.5225, 0.05, 0.38, 0.015])
    a = fig.colorbar(plot2, cax=cbar_ax1, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin, caxismax, 9))
    a.set_label(label="Change in amplitude of divergence of u'v' [m/s/s x 1e6]",
                size=20, weight='bold')
    cl = plt.getp(cbar_ax1, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cols = ['Latitude', 'Magnitude']
    rows = ['Winter response', 'Summer response']
    pad = 10  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline',
                    fontweight='bold', fontsize=24)
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90,
                    fontweight='bold', fontsize=24)
    plt.show()


def zonalvTsens(data, p, sigma, lat, control, heat, csigma, clat):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    plt.suptitle("v'T' sensitivity", size=30, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscaleint(5)
    caxismin1, caxismax1, ctr1 = colourscaleint(10)
    cctrs = [-25, -20, -15, -10, -5, 5, 10]

    plot = axs[0, 0].contourf(meshlat, meshsigma, data[0, 1, :], ctr1,
                              extend='both',
                              cmap=newcmap, vmin=caxismin1, vmax=caxismax1)
    plot1 = axs[0, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot = axs[0, 1].contourf(meshlat, meshsigma, data[0, 0, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plot1 = axs[0, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot = axs[1, 0].contourf(meshlat, meshsigma, data[1, 1, :], ctr1,
                              extend='both', cmap=newcmap, vmin=caxismin1,
                              vmax=caxismax1)
    plot1 = axs[1, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot2 = axs[1, 1].contourf(meshlat, meshsigma, data[1, 0, :], ctrs,
                               extend='both', cmap=newcmap, vmin=caxismin,
                               vmax=caxismax)
    plot1 = axs[1, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    for i in range(2):
        for j in range(2):
            axs[i, j].invert_yaxis()
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(),
                                      fontweight='bold', fontsize=20)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(),
                                      fontweight='bold', fontsize=20)
            meshsigma2 = np.copy(meshsigma)
            meshsigma1 = np.copy(meshsigma)
            for a in range(9):
                for b in range(34):
                    if p[i, 1-j, a, b] <= 0.005:
                        meshsigma2[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma2, marker='o', color='gray',
                              s=0.5)
            for a in range(9):
                for b in range(34):
                    if p[i, 1-j, a, b] > 0.005:
                        meshsigma1[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma1, marker='o', color='black',
                              s=5)
    plt.subplots_adjust(hspace=0.079, wspace=.05)
    cbar_ax = fig.add_axes([0.125, 0.05, 0.38, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin1, caxismax1, 9))
    b.set_label(label="Poleward shift in v'T' [deg]", size=20, weight='bold')
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cbar_ax1 = fig.add_axes([0.5225, 0.05, 0.38, 0.015])
    a = fig.colorbar(plot2, cax=cbar_ax1, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin, caxismax, 9))
    a.set_label(label="Amplitude change of v'T' [K m/s]", size=20,
                weight='bold')
    cl = plt.getp(cbar_ax1, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cols = ['Latitude', 'Magnitude']
    rows = ['Winter response', 'Summer response']
    pad = 10  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline',
                    fontweight='bold', fontsize=24)
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90,
                    fontweight='bold', fontsize=24)
    plt.show()


def zonaluspread(data, sigma, lat, control, csigma, clat):
    '''
    code for calculating std, instead of percentiles
    jetspread_std=np.zeros((2,2,307))
    for i in range(2):
        for j in range(2):
            for k in range(307):
                jetspread_std[i,j,k]=np.std(jetind[k,:,i,j])
                jetspread_std[i,j,k]=jetspread_std[i,j,k]-jetspread_std[i,j,0]
    jetspread_std=jetspread_std[:,:,1:]
    jetspread_std=np.reshape(jetspread_std,(2,2,9,34))
    '''
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    plt.suptitle("Jet spread sensitivity", size=30, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscaleint1(0.5)
    caxismin1, caxismax1, ctr1 = colourscaleint1(5)
    cmi, cma, cctrs = colourscaleint(50)

    plot = axs[0, 0].contourf(meshlat, meshsigma, data[0, 1, :], ctr1,
                              extend='both',
                              cmap=newcmap, vmin=caxismin1, vmax=caxismax1)
    plot1 = axs[0, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.f')

    plot = axs[0, 1].contourf(meshlat, meshsigma, data[0, 0, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plot1 = axs[0, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[0, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    plot = axs[1, 0].contourf(meshlat, meshsigma, data[1, 1, :], ctr1,
                              extend='both', cmap=newcmap, vmin=caxismin1,
                              vmax=caxismax1)
    plot1 = axs[1, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 0].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.f')

    plot2 = axs[1, 1].contourf(meshlat, meshsigma, data[1, 0, :], ctrs,
                               extend='both', cmap=newcmap, vmin=caxismin,
                               vmax=caxismax)
    plot1 = axs[1, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    axs[1, 1].clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                     fmt='%.0f')

    for i in range(2):
        for j in range(2):
            axs[i, j].invert_yaxis()
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(),
                                      fontweight='bold', fontsize=20)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(),
                                      fontweight='bold', fontsize=20)

    plt.subplots_adjust(hspace=0.079, wspace=.05)
    cbar_ax = fig.add_axes([0.125, 0.05, 0.38, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin1, caxismax1, 5))
    b.set_label(label='Increase in standard deviation of jet latitude index [deg]', size=20,
                weight='bold', fontsize=20)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cbar_ax1 = fig.add_axes([0.5225, 0.05, 0.38, 0.015])
    a = fig.colorbar(plot2, cax=cbar_ax1, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin, caxismax, 5))
    a.set_label(label='Increase in standard deviation of jet speed index [m/s]', size=20,
                weight='bold', fontsize=20)
    cl = plt.getp(cbar_ax1, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cols = ['Latitude', 'Magnitude']
    rows = ['Winter response', 'Summer response']
    pad = 10  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline',
                    fontweight='bold', fontsize=24)
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90,
                    fontweight='bold', fontsize=24)
    plt.show()


def zonalall(u, duv, potemp, strm, vT, divvT, temp, z, sigma, lat, heat, run):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    fig, axs = plt.subplots(2, 2, sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Run " + str(run), size=30, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    lapse = lapserate(temp[run, :, :], z[run, :, :], sigma, lat)
    lapse[22:, :] = 99
    lapse0 = lapserate(temp[0, :, :], z[0, :, :], sigma, lat)
    lapse0[22:, :] = 99

    ctrsa = np.arange(-65, 0, 5)
    ctrsb = np.arange(5, 70, 5)
    ctrs = np.concatenate((ctrsa, ctrsb))
    ctrsa1 = np.arange(-70, 0, 7)
    ctrsb1 = np.arange(7, 70, 7)
    ctrs2 = np.concatenate((ctrsa1, ctrsb1))
    plot2 = axs[0, 0].contourf(meshlat, meshsigma, u[run, :], ctrs,
                               cmap=newcmap, vmin=-65, vmax=65, extend='both')
    m = axs[0, 0].get_figure().colorbar(plot2, ax=axs[0, 0],
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    m.set_label(label="Zonal wind [m/s]\n(Contours show divergence of u'v' [m/s/s x 1e6])", weight='bold', fontsize=20)
    q = m.ax.get_xticklabels()
    m.ax.set_xticklabels(q, fontweight='bold', fontsize=20)
    axs[0, 0].contour(meshlat, meshsigma, duv[0, :]*1e6, ctrs2,
                      colors='g', linewidths=1.5)
    plot1 = axs[0, 0].contour(meshlat, meshsigma, duv[run, :]*1e6, ctrs2,
                              colors='k', linewidths=1.5)
    plt.clabel(plot1, inline=True, inline_spacing=0, fontsize=10, fmt='%.0f')

    ctrsa = np.arange(-2, 0, .2)
    ctrsb = np.arange(0, 2.2, .2)
    ctrs = np.concatenate((ctrsa, ctrsb))
    #ctrs2 = [270, 280, 290, 300, 310, 320, 330, 340, 350, 400, 500]
    ctrs2 = np.linspace(270, 350, 17)
    plot3 = axs[1, 0].contourf(meshlat, meshsigma, strm[run, :]/1e11, ctrs,
                               cmap=newcmap, vmin=-2, vmax=2, extend='both')
    n = axs[0, 0].get_figure().colorbar(plot3, ax=axs[1, 0],
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    n.set_label(label="Streamfunction [kg/s x 1e-11]\n(Contours show potential temperature [K])", weight='bold', fontsize=20)
    r = n.ax.get_xticklabels()
    n.ax.set_xticklabels(r, fontweight='bold', fontsize=20)
    axs[1, 0].contour(meshlat, meshsigma, potemp[0, :], ctrs2,
                      colors='g', linewidths=1.5)
    plot1 = axs[1, 0].contour(meshlat, meshsigma, potemp[run, :], ctrs2,
                              colors='k', linewidths=1.5)
    axs[1, 0].clabel(plot1, inline=True, inline_spacing=-10, fontsize=10,
                     fmt='%.1f')

    cmin, cmax, ctrs = colourscale(u[run, :]-u[0, :])
    ctrsc1 = np.arange(-30, 0, 3)
    ctrsd1 = np.arange(3, 30, 3)
    ctrs2 = np.concatenate((ctrsc1, ctrsd1))
    plot4 = axs[0, 1].contourf(meshlat, meshsigma, u[run, :]-u[0, :], ctrs,
                               cmap=newcmap, vmin=cmin, vmax=cmax,
                               extend='both')
    o = axs[0, 1].get_figure().colorbar(plot4, ax=axs[0, 1],
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    o.set_label(label="Zonal wind change [m/s]\n(Contours show v'T' [K m/s])",
                weight='bold', fontsize=20)
    s = o.ax.get_xticklabels()
    o.ax.set_xticklabels(s, fontweight='bold', fontsize=20)
    axs[0, 1].contour(meshlat, meshsigma, vT[0, :], ctrs2, colors='g',
                      linewidths=1.5)
    plot1 = axs[0, 1].contour(meshlat, meshsigma, vT[run, :], ctrs2,
                              colors='k', linewidths=1.5)
    axs[0, 1].clabel(plot1, inline=True, inline_spacing=0, fontsize=10,
                     fmt='%.0f')

    cmin, cmax, ctrs = colourscale((strm[run, :]-strm[0, :])/1e10)
    ctrsc = np.arange(-21, 0, 3)
    ctrsd = np.arange(3, 24, 3)
    ctrs2 = np.concatenate((ctrsc, ctrsd))
    plot5 = axs[1, 1].contourf(meshlat, meshsigma,
                               (strm[run, :]-strm[0, :])/1e10, ctrs,
                               cmap=newcmap, vmin=cmin, vmax=cmax,
                               extend='both')
    p = axs[1, 1].get_figure().colorbar(plot5, ax=axs[1, 1],
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    t = p.ax.get_xticklabels()
    p.ax.set_xticklabels(t, fontweight='bold', fontsize=20)
    p.set_label(label="Streamfunction change [kg/s x 1e-10]\n(Contours show divergence of v'T' [K/s x 1e6])", weight='bold', fontsize=20)    
    axs[1, 1].contour(meshlat, meshsigma, divvT[0, :]*1e6, ctrs2,
                      colors='g', linewidths=1.5)
    plot1 = axs[1, 1].contour(meshlat, meshsigma, divvT[run, :]*1e6, ctrs2,
                              colors='k', linewidths=1.5)
    axs[1, 1].clabel(plot1, inline=True, inline_spacing=0, fontsize=10,
                     fmt='%.0f')
    for i in range(2):
        for j in range(2):
            axs[i, j].contour(meshlat[1:, :], meshsigma[1:, :], lapse0, [2])
            axs[i, j].contour(meshlat[1:, :], meshsigma[1:, :], lapse, [2],
                              linewidths=2)
            axs[i, j].invert_yaxis()
            if run != 0:
                sigmax = np.argmax(heat[run, :, :], axis=0)[0]
                latmax = np.argmax(heat[run, :, :], axis=1)[0]
                axs[i, j].scatter(lat[latmax], sigma[sigmax], marker='o',
                                  color='deeppink', linewidth=2, s=50)
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(),
                                      fontweight='bold', fontsize=20)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(),
                                      fontweight='bold', fontsize=20)
    plt.subplots_adjust(hspace=0.15, wspace=.05)
    plt.show()


def comparison():
    def lapserate(t, z, sigma, lat):
        """
        Produces plot of lapse rate of T data

        Parameters
        ----------
        t: array
            temperature data field
        z: array
            geopotential height of field
        sigma: array
            sigma levels of field
        lat: array
            lat of field
        """
        import numpy as np
        dT = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
        dz = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
        for i in range(np.ma.size(sigma, axis=0)-1):
            dT[i, :] = t[i+1, :] - t[i, :]
        for i in range(np.ma.size(sigma, axis=0)-1):
            dz[i, :] = z[i+1, :] - z[i, :]
        lapse = -1000 * dT[0:-1] / dz[0:-1]
        # zonalplot(lapse, sigma[0:-1], lat, 'Lapse rate')
        return lapse

    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from adjustText import adjust_text
    meshlat000, meshsigma = np.meshgrid(lat, sigma)
    meshlat, meshlevel = np.meshgrid(lat000,level)
    lapse = lapserate(airNA, hgtNA, level, lat000)
    lapsegcm = lapserate(temp, z000[0, :], sigma, lat)
    lapse[:4, :] = 99
    lapsegcm[22:, :] = 99
    fig = plt.figure()
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscaleint1(50)
    plt.subplot(1, 2, 1)
    plot = plt.contourf(meshlat, meshlevel, uNA, ctrs, extend='both',
                        cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plt.xticks(weight='bold', fontsize=20)
    plt.yticks(weight='bold', fontsize=20)
    contour(meshlat[1:, :], meshlevel[1:, :], lapse, [2], linewidths=2)
    plt.axvline(0, color='dimgrey', linestyle='dashed')
    plt.gca().invert_yaxis()
    plt.yscale('linear')
    plt.xlim([-90, 90])
    plt.xticks(np.arange(-75, 90, 15))
    plt.title('Reanalysis', y=1.02, fontsize=30, weight='bold')
    plt.subplot(1, 2, 2)
    plot = plt.contourf(meshlat000, meshsigma, u, ctrs, extend='both',
                        cmap=newcmap, vmin=caxismin, vmax=caxismax)
    CS = plt.contour(lat, sigma, heat[124, :],
                     [0.5, 1.0, 1.5, 1.999], colors='k', linewidths=3)
    plt.clabel(CS, inline=True, inline_spacing=-5, fontsize=11,
               fmt='%.1f', colors='k', weight='bold')
    for i in range(37):
        plt.axhline(y=sigma[i], xmin=0.025, xmax=0.075, color='k', lw=2)
    plt.xticks(weight='bold', size=20)
    contour(meshlat000[1:, :], meshsigma[1:, :], lapsegcm, [2], linewidths=2)
    
    sigmax = np.argmax(heat[62, :, :], axis=0)[0]
    latmax = np.argmax(heat[62, :, :], axis=1)[0]
    plt.scatter(lat[latmax], sigma[sigmax], marker='o', 
                    color='deeppink', linewidth=2, s=50)
    sigmax = np.argmax(heat[160, :, :], axis=0)[0]
    latmax = np.argmax(heat[160, :, :], axis=1)[0]
    plt.scatter(lat[latmax], sigma[sigmax], marker='o', 
                    color='deeppink', linewidth=2, s=50)
    sigmax = np.argmax(heat[175, :, :], axis=0)[0]
    latmax = np.argmax(heat[175, :, :], axis=1)[0]
    plt.scatter(lat[latmax], sigma[sigmax], marker='o', 
                    color='deeppink', linewidth=2, s=50)
    sigmax = np.argmax(heat[215, :, :], axis=0)[0]
    latmax = np.argmax(heat[215, :, :], axis=1)[0]
    plt.scatter(lat[latmax], sigma[sigmax], marker='o', 
                    color='deeppink', linewidth=2, s=50)
    labels1 = ['1', '2', '3', '4']
    texts = []
    for x, y, s in zip([54.4, 32.1, -65.6, -32.1],
                       [0.9,0.6, 0.5, 0.4],
                       labels1):
        texts.append(plt.text(x, y, s, bbox={'pad': 0, 'alpha': 0}, size=25,
                     weight='bold'))
    adjust_text(texts, expand_points=(0, 0),
                force_points=0.01, arrowprops=dict(arrowstyle="->", color='k', lw=0),
                bbox={'pad': 5, 'alpha': 0}, size=10)
    plt.gca().invert_yaxis()
    plt.yscale('linear')
    plt.xlim([-90, 90])
    plt.ylim([1, .01])
    plt.xticks(np.arange(-75, 90, 15))
    plt.tick_params(axis='y', labelleft='off')
    plt.title('Model', y=1.02, fontsize=30, weight='bold')
    cbar_ax = fig.add_axes([0.325, 0.05, 0.38, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='horizontal')
    b.set_label(label='Zonal wind [m/s]', size=20,
                weight='bold', fontsize=20)
    plt.subplots_adjust(hspace=0.05, wspace=.05)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    plt.show()
    plt.figure()
    plt.title('Surface temperature', y=1.02, fontsize=30, weight='bold')
    plt.plot(latgauss, surfNA-273.15, color='red', linewidth=2)
    plt.plot(lat, temp[-1, :]-273.15, color='blue', linewidth=2)
    plt.axvline(0, color='dimgrey', linestyle='dashed')
    plt.xlabel('Latitude [deg]', fontweight='bold', fontsize=20)
    plt.ylabel('Temperature [Celsius]', fontweight='bold', fontsize=20)
    plt.yticks(weight='bold', size=20)
    plt.xticks(weight='bold', size=20)
    plt.xlim([-90, 90])
    plt.xticks(np.arange(-75, 90, 15))
    plt.ylim([-35, 30])


def zonalhadsens(data1, p, sigma, lat, control, csigma, clat):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    data = np.copy(data1)
    fig, axs = plt.subplots(2, 2, sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    plt.suptitle("Hadley Cell Sensitivity", size=36, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    data[0, :] = data[0, :] / 1e10
    caxismin, caxismax, ctrs = colourscaleint(3)
    caxismin1, caxismax1, ctr1 = colourscaleint(5)
    control = control / 1e11
    cctrs = [-1.3, -1.2, -1.1, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
             .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2]

    plot = axs[0, 0].contourf(meshlat, meshsigma, data[1, 0, :], ctr1,
                              extend='both',
                              cmap=newcmap, vmin=caxismin1, vmax=caxismax1)
    plot1 = axs[0, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)
    #  axs[0, 0].clabel(plot1, [-1,-.5,.5,1], inline=True, inline_spacing=-3,
    #                   fontsize=10,
    #                   fmt='%.1f')

    plot = axs[0, 1].contourf(meshlat, meshsigma, data[0, 0, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plot1 = axs[0, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)

    plot = axs[1, 0].contourf(meshlat, meshsigma, data[1, 1, :], ctr1,
                              extend='both', cmap=newcmap, vmin=caxismin1,
                              vmax=caxismax1)
    plot1 = axs[1, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)

    plot2 = axs[1, 1].contourf(meshlat, meshsigma, data[0, 1, :], ctrs,
                               extend='both', cmap=newcmap, vmin=caxismin,
                               vmax=caxismax)
    plot1 = axs[1, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)

    for i in range(2):
        for j in range(2):
            axs[i, j].invert_yaxis()
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0.1)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.2, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(),
                                      fontweight='bold', fontsize=20)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(),
                                      fontweight='bold', fontsize=20)
            meshsigma2 = np.copy(meshsigma)
            meshsigma1 = np.copy(meshsigma)
            # axs[i, j].contour(meshlat, meshsigma, p[i, 1-j, :], [0.005],
            # linewidths=2)
            for a in range(9):
                for b in range(34):
                    if p[1-j, i, a, b] <= 0.005:
                        meshsigma2[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma2, marker='o', color='gray',
                              s=0.5)
            for a in range(9):
                for b in range(34):
                    if p[1-j, i, a, b] > 0.005:
                        meshsigma1[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma1, marker='o', color='black',
                              s=5)
    plt.subplots_adjust(hspace=0.079, wspace=.05)
    cbar_ax = fig.add_axes([0.125, 0.05, 0.38, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin1, caxismax1, 9))
    b.set_label(label='Poleward Hadley Cell extent shift [deg]', size=20,
                weight='bold')
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cbar_ax1 = fig.add_axes([0.5225, 0.05, 0.38, 0.015])
    a = fig.colorbar(plot2, cax=cbar_ax1, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin, caxismax, 9))
    a.set_label(label='Change in strength of Hadley Cell [kg/s x 1e10]', size=20, weight='bold')
    cl = plt.getp(cbar_ax1, 'xmajorticklabels')
    plt.setp(cl, fontweight='bold', fontsize=20)
    cols = ['Latitude', 'Magnitude']
    rows = ['Winter response', 'Summer response']
    pad = 10  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline',
                    fontweight='bold', fontsize=24)
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90,
                    fontweight='bold', fontsize=24)

    plt.show()


def zonalHad_case(strm, sigma, lat, heat, run):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    # plt.suptitle("Run " + str(run), size=30, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[242:270, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    data = (strm[run, :] - strm[0, :]) / 1e10
    caxismin, caxismax, ctrs = colourscaleint(2)
    datacontrol = strm[0, :] / 1e11
    caxismin1, caxismax1, ctr1 = colourscaleint(2)
    cctrs = [-1.3, -1.2, -1.1, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
             .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2]

    plot2 = plt.contourf(meshlat, meshsigma, data, ctrs,
                         cmap=newcmap, vmin=caxismin, vmax=caxismax,
                         extend='both')
    m = plt.colorbar(plot2, orientation='horizontal', aspect=50,
                     format='%.1f', spacing='proportional')
    m.set_label(label="Streamfunction change [kg/s x 1e10]\n(Contours show streamfunction and control (black and green) [kg/s x 1e11])", weight='bold', fontsize=20)
    q = m.ax.get_xticklabels()
    m.ax.set_xticklabels(q, fontweight='bold', fontsize=20)
    plt.contour(meshlat, meshsigma, datacontrol, cctrs,
                colors='g', linewidths=1.5)
    plt.contour(meshlat, meshsigma, strm[run, :]/1e11, cctrs,
                colors='k', linewidths=1.5)
    ax = plt.gca()
    ax.invert_yaxis()
    if run != 0:
        sigmax = np.argmax(heat[run, :, :], axis=0)[0]
        latmax = np.argmax(heat[run, :, :], axis=1)[0]
        plt.scatter(lat[latmax], sigma[sigmax], marker='o',
                    color='deeppink', linewidth=2, s=50)
    ax.set_xlim([lat[0], lat[-1]])
    ax.set_ylim(1, 0)
    ax.set_xticks(np.arange(-75, 90, 15))
    ax.set_yticks(np.arange(0.0, 1.2, 0.2))
    ax.set_xticklabels(ax.get_xticks(), fontweight='bold', fontsize=20)
    ax.set_yticklabels(ax.get_yticks(), fontweight='bold', fontsize=20)
    plt.show()


def hadscatter(strklats, strklats_actual):
    import matplotlib.pyplot as plt
    plt.figure()
    a = np.reshape(strklats[0, 1:], (9, 34))
    b = np.reshape(a[:, :21], (189))
    c = np.reshape(a[:, 21:], (117))
    d = np.reshape(strklats_actual[0, 1:], (9, 34))
    e = np.reshape(d[:, :21], (189))
    f = np.reshape(d[:, 21:], (117))
    plt.scatter(b, e, color='b')
    plt.scatter(c, f, color='r')
    plt.scatter(strklats[0, 0], strklats_actual[0, 0], color='g')
    plt.plot([-100, 100], [-100, 100], c='k', linewidth=2, alpha=0.5)
    mx = np.max([np.max(strklats[0, :]), np.max(strklats_actual[0, :])])
    mn = np.min([np.min(strklats[0, :]), np.min(strklats_actual[0, :])])
    axes = plt.gca()
    axes.set_xlim([mn, mx])
    axes.set_ylim([mn, mx])
    plt.title('Scatter of actual storm track latitude (y-axis) against predicted latitudes (x-axis)')
    plt.show()
    plt.figure()
    a = np.reshape(strklats[1, 1:], (9, 34))
    b = np.reshape(a[:, :21], (189))
    c = np.reshape(a[:, 21:], (117))
    d = np.reshape(strklats_actual[1, 1:], (9, 34))
    e = np.reshape(d[:, :21], (189))
    f = np.reshape(d[:, 21:], (117))
    plt.scatter(b, e, color='b')
    plt.scatter(c, f, color='r')
    plt.scatter(strklats[1, 0], strklats_actual[1, 0], color='g')
    plt.plot([-100, 100], [-100, 100], c='k', linewidth=2, alpha=0.5)
    mx = np.max([np.max(strklats[1, :]), np.max(strklats_actual[1, :])])
    mn = np.min([np.min(strklats[1, :]), np.min(strklats_actual[1, :])])
    axes = plt.gca()
    axes.set_xlim([mn, mx])
    axes.set_ylim([mn, mx])
    plt.title('Scatter of actual storm track latitude (y-axis) against predicted latitudes (x-axis)')
