#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 11:03:09 2018

@author: bakerh
"""
import numpy as np


def nao_plot(eofhad, eofncep, eofncar, eofera, winter='yes'):
    from mpl_toolkits.basemap import Basemap, addcyclic
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    eofs = [eofhad, eofncep, eofncar, eofera]
    label = ['a', 'b', 'c', 'd']
    titles = ['HadAM3P', 'NCEP2', 'NOAA20C', 'ERA20C']
    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)
    if winter == 'yes':
        colormax = 5
        colormin = -colormax
        ws = 0
        var = ['38.3', '48.7', '41.7', '41.8']
    else:
        colormax = 2.5
        colormin = -colormax
        ws = 1
        var = ['28.2', '43.1', '33.8', '33.2']

    fig, axs = plt.subplots(2, 2, facecolor='w', edgecolor='k', linewidth=2)
    for i in range(4):
        ax1 = axs[int(np.floor(i/2)), np.remainder(i, 2)]
        pdata = eofs[i][ws]
        pdata, lon1 = addcyclic(pdata, lon)
        meshlon, meshlat = np.meshgrid(lon1, lat)
        m = Basemap(width=10000000, height=7000000,
                    resolution='c', projection='aea',
                    lat_1=40., lat_2=60, lon_0=-25, lat_0=50, ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        ctrs = np.linspace(colormin, colormax, 17)
        plot = m.contourf(x, y, pdata, ctrs,
                          cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                          extend='both')
        parallels = m.drawparallels(np.arange(-90., 75., 15.))
        meridians = m.drawmeridians(np.arange(-180., 181., 30))
        m.drawparallels(parallels, labels=[True, False, False, False])
        m.drawmeridians(meridians, labels=[False, False, False, True])
        ax1.annotate(label[i], xy=(-.04, .95), fontname='Arial', fontsize=20,
                     fontweight='bold', xycoords='axes fraction')
        ax1.annotate('Var = ' + var[i] + '%', xy=(.78, 1.02), fontname='Arial',
                     fontsize=16, xycoords='axes fraction')
        ax1.set_title(titles[i], fontname='Arial', fontsize=16)

    cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.01])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.2f')
    b.set_label(label='NAO (hPa per SD)', size=16,
                fontsize=16, fontname='Arial', labelpad=0)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    plt.subplots_adjust(hspace=0.2, wspace=0.08, top=.95, bottom=0.15,
                        left=0.03, right=.97)
    plt.show()


def mapplotssteof(eofs, mx=1, mask='yes', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    long = np.concatenate((lon, [360]))
    titles = ['EOF1', 'EOF2', 'EOF3', 'EOF4']
    var = ['60.2', '10.2', '7.6', '4.3']
    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(4):
        ax1 = axs[int(np.floor(a/2)), np.remainder(a, 2)]
        plotdata1 = np.concatenate((eofs[a],
                                    np.expand_dims(eofs[a, :, 0],
                                    axis=1)), axis=1)
        if mask == 'yes':
            lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
            lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                 axis=1)
            plotdata1 = np.ma.masked_array(plotdata1, lsm)
        meshlon, meshlat = np.meshgrid(long, lat)
        m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                    llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        mycmap2 = plt.cm.YlOrRd(np.arange(256))
        mycmap1 = plt.cm.Blues_r(np.arange(256))
        my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
        my_cmap[238:275, :] = 1
        newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                               my_cmap)
        caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
        plot = m.contourf(x, y, plotdata1, ctrs,
                          cmap=newcmap, vmin=caxismin, vmax=caxismax,
                          extend='both')
        ax1.set_ylim(-60, 60)
        m.drawparallels(np.arange(-60., 90., 30.),
                        labels=[True, False, False, True], linewidth=0)
        m.drawmeridians(np.arange(0., 390., 30.),
                        labels=[True, False, False, True], linewidth=0)
        ax1.set_yticks(np.arange(-60., 90., 30.))
        ax1.set_xticks(np.arange(0., 390., 30.))
        ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                        direction='out', length=5, width=2)
        ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
        ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
        ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                        direction='out', length=4, width=1)
        ax1.annotate('Var = ' + var[a] + '%', xy=(.86, 1.02), fontname='Arial',
                     fontsize=16, xycoords='axes fraction')
        ax1.set_title(titles[a], fontname='Arial', fontsize=16)
        poly1 = mpl.patches.Polygon([m(30, -30), m(105, -30),
                                    m(105, 20), m(30, 20)],
                                   linewidth=2, edgecolor='fuchsia', fill=None)
        ax1.add_patch(poly1)
    cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.01])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.2f')
    b.set_label(label='SST (K per SD)', size=16,
                fontsize=16, fontname='Arial', labelpad=0)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    plt.subplots_adjust(hspace=0.2, wspace=0.08, top=.95, bottom=0.15,
                        left=0.03, right=.97)
    plt.show()


def nao_plot_period(eof, eof1, eof2, eof3, winter='no'):
    from mpl_toolkits.basemap import Basemap, addcyclic
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    eofs = [eof, eof1, eof2, eof3]
    label = ['a', 'b', 'c', 'd']
    titles = ['1872-2012', '1872-1939', '1940-1979', '1980-2012']
    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)
    if winter == 'yes':
        colormax = 5
        colormin = -colormax
        ws = 0
        var = ['41.7', '36.0', '42.8', '51.8']
    else:
        colormax = 2.5
        colormin = -colormax
        ws = 1
        var = ['33.8', '36.0', '31.1', '38.0']

    fig, axs = plt.subplots(2, 2, facecolor='w', edgecolor='k', linewidth=2)
    for i in range(4):
        ax1 = axs[int(np.floor(i/2)), np.remainder(i, 2)]
        pdata = eofs[i][ws]
        pdata, lon1 = addcyclic(pdata, lon)
        meshlon, meshlat = np.meshgrid(lon1, lat)
        m = Basemap(width=10000000, height=7000000,
                    resolution='c', projection='aea',
                    lat_1=40., lat_2=60, lon_0=-25, lat_0=50, ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        ctrs = np.linspace(colormin, colormax, 17)
        plot = m.contourf(x, y, pdata, ctrs,
                          cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                          extend='both')
        parallels = m.drawparallels(np.arange(-90., 75., 15.))
        meridians = m.drawmeridians(np.arange(-180., 181., 30))
        m.drawparallels(parallels, labels=[True, False, False, False])
        m.drawmeridians(meridians, labels=[False, False, False, True])
        ax1.annotate(label[i], xy=(-.04, .95), fontname='Arial', fontsize=20,
                     fontweight='bold', xycoords='axes fraction')
        ax1.annotate('Var = ' + var[i] + '%', xy=(.78, 1.02), fontname='Arial',
                     fontsize=16, xycoords='axes fraction')
        ax1.set_title(titles[i], fontname='Arial', fontsize=16)

    cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.01])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.2f')
    b.set_label(label='NAO (hPa per SD)', size=16,
                fontsize=16, fontname='Arial', labelpad=0)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    plt.subplots_adjust(hspace=0.2, wspace=0.08, top=.95, bottom=0.15,
                        left=0.03, right=.97)
    plt.show()


def ea_plot(eofhad, eofncep, eofncar, eofera, winter='yes'):
    from mpl_toolkits.basemap import Basemap, addcyclic
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    eofs = [eofhad, eofncep, eofncar, eofera]
    label = ['a', 'b', 'c', 'd']
    titles = ['HadAM3P', 'NCEP2', 'NOAA20C', 'ERA20C']
    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)
    if winter == 'yes':
        colormax = 5
        colormin = -colormax
        ws = 0
        var = ['16.9', '10.1', '18.4', '18.0']
    else:
        colormax = 2
        colormin = -colormax
        ws = 1
        var = ['14.8', '12.3', '13.5', '16.2']

    fig, axs = plt.subplots(2, 2, facecolor='w', edgecolor='k', linewidth=2)
    for i in range(4):
        ax1 = axs[int(np.floor(i/2)), np.remainder(i, 2)]
        pdata = eofs[i][ws]
        pdata, lon1 = addcyclic(pdata, lon)
        meshlon, meshlat = np.meshgrid(lon1, lat)
        m = Basemap(width=10000000, height=7000000,
                    resolution='c', projection='aea',
                    lat_1=40., lat_2=60, lon_0=-25, lat_0=50, ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        ctrs = np.linspace(colormin, colormax, 17)
        plot = m.contourf(x, y, pdata, ctrs,
                          cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                          extend='both')
        parallels = m.drawparallels(np.arange(-90., 75., 15.))
        meridians = m.drawmeridians(np.arange(-180., 181., 30))
        m.drawparallels(parallels, labels=[True, False, False, False])
        m.drawmeridians(meridians, labels=[False, False, False, True])
        ax1.annotate(label[i], xy=(-.04, .95), fontname='Arial', fontsize=20,
                     fontweight='bold', xycoords='axes fraction')
        ax1.annotate('Var = ' + var[i] + '%', xy=(.78, 1.02), fontname='Arial',
                     fontsize=16, xycoords='axes fraction')
        ax1.set_title(titles[i], fontname='Arial', fontsize=16)

    cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.01])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.2f')
    b.set_label(label='EA (hPa per SD)', size=16,
                fontsize=16, fontname='Arial', labelpad=0)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    plt.subplots_adjust(hspace=0.2, wspace=0.08, top=.95, bottom=0.15,
                        left=0.03, right=.97)
    plt.show()


def mapplot60(plotdata, plotdata2, ji, mask='yes', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(2):
        for b in range(2):
            if b == 0:
                mx = 1.2
            else:
                mx = .8
            plotdata1 = np.concatenate((plotdata[a, 1-b, :],
                                        np.expand_dims(plotdata[a, 1-b, :, 0],
                                        axis=1)), axis=1)
            plotdata3 = np.concatenate((plotdata2[a, 1-b, :],
                                        np.expand_dims(plotdata2[a, 1-b, :, 0],
                                        axis=1)), axis=1)
            if mask == 'yes':
                lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
                lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                     axis=1)
                plotdata1 = np.ma.masked_array(plotdata1, lsm)
                plotdata3 = np.ma.masked_array(plotdata3, lsm)
            meshlon, meshlat = np.meshgrid(long, lat)
            ax1 = axs[a, b]
            m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                        llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
            m.drawcoastlines()
            m.drawmapboundary(linewidth=2)
            x, y = m(meshlon, meshlat)
            mycmap2 = plt.cm.YlOrRd(np.arange(256))
            mycmap1 = plt.cm.Blues_r(np.arange(256))
            my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
            #my_cmap[239:274, :] = 1
            newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                                   my_cmap)
            caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
            plot = m.contourf(x, y, plotdata1, ctrs,
                              cmap=newcmap, vmin=caxismin, vmax=caxismax,
                              extend='both')
            m.plot(lon[160:], ji[a, 160:], color='red')
            m.contour(x, y, plotdata3, ctrs, colors='k')
            ax1.set_ylim(-60, 60)
            m.drawparallels(np.arange(-60, 90, 30),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0, 390, 30),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-60., 90., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if a == 0 and b == 0:
                ax1.set_title('Winter Latitude', fontsize=16, y=1.03)
            if a == 1 and b == 0:
                ax1.set_title('Summer Latitude', fontsize=16, y=1.03)
                c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                                 spacing='proportional', aspect=50)
                c.set_label(label='Poleward jet latitude shift ($^\circ$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
            if a == 0 and b == 1:
                ax1.set_title('Winter Speed', fontsize=16, y=1.03)
            if a == 1 and b == 1:
                ax1.set_title('Summer Speed', fontsize=16, y=1.03)
                c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                                 spacing='proportional', aspect=50)
                c.set_label(label='Jet speed increase (ms$\mathregular{^{-1}}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    plt.subplots_adjust(hspace=0, wspace=0.1, top=.95, bottom=0.05, left=.05,
                        right=.95)

    plt.show()


def mapplot60_nao(plotdata_sig, plotdata, mx=1, mask='yes', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(2, 1, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(2):
        ax1 = axs[a]
        plotdata1 = np.concatenate((plotdata_sig[a],
                                    np.expand_dims(plotdata[a, :, 0],
                                    axis=1)), axis=1)
        plotdata2 = np.concatenate((plotdata[a],
                                    np.expand_dims(plotdata[a, :, 0],
                                    axis=1)), axis=1)
        if mask == 'yes':
            lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
            lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                 axis=1)
            plotdata1 = np.ma.masked_array(plotdata1, lsm)
            plotdata2 = np.ma.masked_array(plotdata2, lsm)
        meshlon, meshlat = np.meshgrid(long, lat)
        m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                    llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        mycmap2 = plt.cm.YlOrRd(np.arange(256))
        mycmap1 = plt.cm.Blues_r(np.arange(256))
        my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
        # my_cmap[238:275, :] = 1
        newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                               my_cmap)
        caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
        plot = m.contourf(x, y, plotdata1, ctrs,
                          cmap=newcmap, vmin=caxismin, vmax=caxismax,
                          extend='both')
        m.contour(x, y, plotdata2, ctrs, colors='k')
        ax1.set_ylim(-60, 60)
        m.drawparallels(np.arange(-60., 90., 30.),
                        labels=[True, False, False, True], linewidth=0)
        m.drawmeridians(np.arange(0., 390., 30.),
                        labels=[True, False, False, True], linewidth=0)
        ax1.set_yticks(np.arange(-60., 90., 30.))
        ax1.set_xticks(np.arange(0., 390., 30.))
        ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                        direction='out', length=5, width=2)
        ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
        ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
        ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                        direction='out', length=4, width=1)
        if a == 0:
            ax1.set_ylabel('1872-1939', fontsize=16, labelpad=20)
        else:
            ax1.set_ylabel('1940-1979', fontsize=16, labelpad=20)
    c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                     spacing='proportional', aspect=50)
    c.set_label(label='NAO sensitivity ([$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    #c.set_label(label='Circulation index sensitivity (ms$^{-1}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    plt.subplots_adjust(hspace=0, wspace=0.1, top=.95, bottom=0.05, left=.05,
                        right=.95)
    plt.show()


def mapplot60all(gtojsig, gtoj, gtonaosig, gtonao, ji, mask='yes', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(3, 2, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(3):
        for b in range(2):
            if a == 0:
                mx = 1  
                ax1 = axs[a, b]
                plotdata1 = np.concatenate((gtonaosig[b],
                                            np.expand_dims(gtonaosig[b, :, 0],
                                            axis=1)), axis=1)
                plotdata2 = np.concatenate((gtonao[b],
                                            np.expand_dims(gtonao[b, :, 0],
                                            axis=1)), axis=1)
                if mask == 'yes':
                    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
                    lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                         axis=1)
                    plotdata1 = np.ma.masked_array(plotdata1, lsm)
                    plotdata2 = np.ma.masked_array(plotdata2, lsm)
                meshlon, meshlat = np.meshgrid(long, lat)
                m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                            llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
                m.drawcoastlines()
                m.drawmapboundary(linewidth=2)
                x, y = m(meshlon, meshlat)
                mycmap2 = plt.cm.YlOrRd(np.arange(256))
                mycmap1 = plt.cm.Blues_r(np.arange(256))
                my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
                # my_cmap[238:275, :] = 1
                newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                                       my_cmap)
                caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
                plot = m.contourf(x, y, plotdata1, ctrs,
                                  cmap=newcmap, vmin=caxismin, vmax=caxismax,
                                  extend='both')
                m.contour(x, y, plotdata2, ctrs, colors='k')
                m.plot(lon[160:], ji[b, 160:], color='red')
                ax1.set_ylim(-60, 60)
                m.drawparallels(np.arange(-60., 90., 30.),
                                labels=[True, False, False, True], linewidth=0)
                m.drawmeridians(np.arange(0., 390., 30.),
                                labels=[True, False, False, True], linewidth=0)
                ax1.set_yticks(np.arange(-60., 90., 30.))
                ax1.set_xticks(np.arange(0., 390., 30.))
                ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                                direction='out', length=5, width=2)
                ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
                ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
                ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                                direction='out', length=4, width=1)
                if b == 1:
                    ax1.set_title('Summer NAO', fontsize=16, y=1.03)
                    cbar_ax = fig.add_axes([0.925, 0.6842, 0.005, 0.235])
                    d = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', orientation='vertical', extend='max')
                    d.set_label(label='NAO sensitivity\n([$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16, fontsize=16)
                else:
                    ax1.set_title('Winter NAO', fontsize=16, y=1.03)
            elif a == 1:
                mx = 1.2
            else:
                mx = .8
            if a != 0:
                plotdata1 = np.concatenate((gtojsig[b, 2-a, :],
                                            np.expand_dims(gtojsig[b, 2-a, :, 0],
                                            axis=1)), axis=1)
                plotdata3 = np.concatenate((gtoj[b, 2-a, :],
                                            np.expand_dims(gtoj[b, 2-a, :, 0],
                                            axis=1)), axis=1)
                if mask == 'yes':
                    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
                    lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                         axis=1)
                    plotdata1 = np.ma.masked_array(plotdata1, lsm)
                    plotdata3 = np.ma.masked_array(plotdata3, lsm)
                meshlon, meshlat = np.meshgrid(long, lat)
                ax1 = axs[a, b]
                m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                            llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
                m.drawcoastlines()
                m.drawmapboundary(linewidth=2)
                x, y = m(meshlon, meshlat)
                mycmap2 = plt.cm.YlOrRd(np.arange(256))
                mycmap1 = plt.cm.Blues_r(np.arange(256))
                my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
                #my_cmap[239:274, :] = 1
                newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                                       my_cmap)
                caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
                plot = m.contourf(x, y, plotdata1, ctrs,
                                  cmap=newcmap, vmin=caxismin, vmax=caxismax,
                                  extend='both')
                m.plot(lon[160:], ji[b, 160:], color='red')
                m.contour(x, y, plotdata3, ctrs, colors='k')
                ax1.set_ylim(-60, 60)
                m.drawparallels(np.arange(-60, 90, 30),
                                labels=[True, False, False, True], linewidth=0)
                m.drawmeridians(np.arange(0, 390, 30),
                                labels=[True, False, False, True], linewidth=0)
                ax1.set_yticks(np.arange(-60., 90., 30.))
                ax1.set_xticks(np.arange(0., 390., 30.))
                ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                                direction='out', length=5, width=2)
                ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
                ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
                ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                                direction='out', length=4, width=1)
                if a == 1 and b == 0:
                    ax1.set_title('Winter Jet Latitude', fontsize=16, y=1.03)
                if a == 1 and b == 1:
                    ax1.set_title('Summer Jet Latitude', fontsize=16, y=1.03)
                    cbar_ax = fig.add_axes([0.925, 0.3825, 0.005, 0.235])
                    d = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max')
                    d.set_label(label='Poleward jet latitude shift\n($^\circ$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16, fontsize=16)
                if a == 2 and b == 0:
                    ax1.set_title('Winter Jet Speed', fontsize=16, y=1.03)
                if a == 2 and b == 1:
                    ax1.set_title('Summer Jet Speed', fontsize=16, y=1.03)
                    cbar_ax = fig.add_axes([0.925, 0.0805, 0.005, 0.235])
                    d = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max')
                    d.set_label(label='Jet speed increase\n(ms$\mathregular{^{-1}}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16, fontsize=16)
    plt.subplots_adjust(hspace=-.2, wspace=0.1, top=.99, bottom=0.01, left=.03,
                        right=.9)

    plt.show()


def mapplot90_wind(plotdata_sig, mx=2, mask='no', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    #plotdata_sig[:, 0, :] = plotdata_sig[:, 1, :]
    #plotdata[:, 0, :] = plotdata[:, 1, :]
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/vwnd.mon.mean.nc','lon')
    lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/vwnd.mon.mean.nc','lat')
    long = np.concatenate((lon, [360.9375]))
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(2, 1, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(2):
        ax1 = axs[a]
        plotdata1 = np.concatenate((plotdata_sig[a],
                                    np.expand_dims(plotdata_sig[a, :, 0],
                                    axis=1)), axis=1)
        #plotdata2 = np.concatenate((plotdata[a],
         #                           np.expand_dims(plotdata[a, :, 0],
          #                          axis=1)), axis=1)
        if mask == 'yes':
            lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
            lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                 axis=1)
            plotdata1 = np.ma.masked_array(plotdata1, lsm)
            #plotdata2 = np.ma.masked_array(plotdata2, lsm)
        meshlon, meshlat = np.meshgrid(long, lat)
        m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        mycmap2 = plt.cm.YlOrRd(np.arange(256))
        mycmap1 = plt.cm.Blues_r(np.arange(256))
        my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
        my_cmap[238:275, :] = 1
        newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                               my_cmap)
        caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
        plot = m.contourf(x, y, plotdata1, ctrs,
                          cmap=newcmap, vmin=caxismin, vmax=caxismax,
                          extend='both')
        #m.contour(x, y, plotdata2, ctrs, colors='k')
        #m.quiver(x[::6,::6], y[::6,::6], plotdata2[::6,::6], plotdata1[::6,::6])
        ax1.set_ylim(-90, 90)
        m.drawparallels(np.arange(-60., 90., 30.),
                        labels=[True, False, False, True], linewidth=0)
        m.drawmeridians(np.arange(0., 390., 30.),
                        labels=[True, False, False, True], linewidth=0)
        ax1.set_yticks(np.arange(-90., 120., 30.))
        ax1.set_xticks(np.arange(0., 390., 30.))
        ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                        direction='out', length=5, width=2)
        ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
        ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
        ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                        direction='out', length=4, width=1)
        polyi = mpl.patches.Polygon([m(45, -20), m(90, -20),
                                    m(90, 00), m(45, 00)],
                                   linewidth=2, edgecolor='fuchsia', fill=None)
        polyp = mpl.patches.Polygon([m(150, -20), m(195, -20),
                                    m(195, 00), m(150, 00)],
                                   linewidth=2, edgecolor='fuchsia', fill=None)
        #if a == 0:
        #    ax1.add_patch(polyi)
        #else:
        #    ax1.add_patch(polyp)
        if a == 0:
            ax1.set_ylabel('Winter', fontsize=16, labelpad=20)
        else:
            ax1.set_ylabel('Summer', fontsize=16, labelpad=20)
    #c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                    # spacing='proportional', aspect=50)
    #c.set_label(label='v200 sensitivity (ms$^{-1}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    cbar_ax = fig.add_axes([0.26, 0.06, 0.48, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='v200 NAO (ms$^{-1}$ s.d.$^{-1}$)', size=16)
    #b.set_label(label='v200 (ms$^{-1}$ s.d.$^{-1}$)', size=16)
    '''
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16) 
    '''
    plt.subplots_adjust(hspace=0.15, wspace=0.1, top=.99, bottom=0.125, left=.05,
                        right=.95)
    plt.show()


def mapplot90_mslp(plotdata_sig, plotdata, mx=1, mask='no', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(2, 1, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(2):
        ax1 = axs[a]
        plotdata1 = np.concatenate((plotdata_sig[a],
                                    np.expand_dims(plotdata[a, :, 0],
                                    axis=1)), axis=1)
        plotdata2 = np.concatenate((plotdata[a],
                                    np.expand_dims(plotdata[a, :, 0],
                                    axis=1)), axis=1)
        if mask == 'yes':
            lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
            lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                 axis=1)
            plotdata1 = np.ma.masked_array(plotdata1, lsm)
            plotdata2 = np.ma.masked_array(plotdata2, lsm)
        meshlon, meshlat = np.meshgrid(long, lat)
        m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        mycmap2 = plt.cm.YlOrRd(np.arange(256))
        mycmap1 = plt.cm.Blues_r(np.arange(256))
        my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
        # my_cmap[238:275, :] = 1
        newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                               my_cmap)
        caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
        ctrs2 = np.linspace(-mx*2, mx*2, 33)
        plot = m.contourf(x, y, plotdata1/100, ctrs,
                          cmap=newcmap, vmin=caxismin, vmax=caxismax,
                          extend='both')
        m.contour(x, y, plotdata2/100, ctrs2, colors='k')
        ax1.set_ylim(-90, 90)
        m.drawparallels(np.arange(-60., 90., 30.),
                        labels=[True, False, False, True], linewidth=0)
        m.drawmeridians(np.arange(0., 390., 30.),
                        labels=[True, False, False, True], linewidth=0)
        ax1.set_yticks(np.arange(-90., 120., 30.))
        ax1.set_xticks(np.arange(0., 390., 30.))
        ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                        direction='out', length=5, width=2)
        ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
        ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
        ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                        direction='out', length=4, width=1)
        polyi = mpl.patches.Polygon([m(45, -20), m(90, -20),
                                    m(90, 00), m(45, 00)],
                                   linewidth=2, edgecolor='fuchsia', fill=None)
        polyp = mpl.patches.Polygon([m(150, -20), m(195, -20),
                                    m(195, 00), m(150, 00)],
                                   linewidth=2, edgecolor='fuchsia', fill=None)
        if a == 0:
            ax1.add_patch(polyi)
        else:
            ax1.add_patch(polyp)
        if a == 0:
            ax1.set_ylabel('Indian', fontsize=16, labelpad=20)
        else:
            ax1.set_ylabel('Pacific', fontsize=16, labelpad=20)
    #c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                    # spacing='proportional', aspect=50)
    #c.set_label(label='v200 sensitivity (ms$^{-1}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    cbar_ax = fig.add_axes([0.26, 0.06, 0.48, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='SLP sensitivity (hPa K$^{-1}$)', size=16)
    '''
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16) 
    '''
    plt.subplots_adjust(hspace=0.15, wspace=0.1, top=.99, bottom=0.125, left=.05,
                        right=.95)
    plt.show()


'''
def mapplot90_windres(eof_had, v200_nao, v_sens):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
    long = np.concatenate((lon, [360.9375]))
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 360, 193)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    eof_2plot = np.concatenate((eof_had[0, :, 96:], eof_had[0, :, :96]), axis=1)
    eof_2plot = np.concatenate((eof_2plot, np.expand_dims(eof_2plot[:, 0],
                                                          axis=1)), axis=1)
    v200_nao = np.concatenate((v200_nao[0], np.expand_dims(v200_nao[0, :, 0],
                               axis=1)), axis=1)
    vi = np.concatenate((v_sens[0], np.expand_dims(v_sens[0, :, 0], axis=1)), axis=1)
    vp = np.concatenate((v_sens[1], np.expand_dims(v_sens[1, :, 0], axis=1)), axis=1)
    proji = np.sum(vi[:, :-1]*v200_nao[:, :-1]*meshlatweight)/np.sum(v200_nao[:, :-1]*v200_nao[:, :-1]*meshlatweight)
    projp = np.sum(vp[:, :-1]*v200_nao[:, :-1]*meshlatweight)/np.sum(v200_nao[:, :-1]*v200_nao[:, :-1]*meshlatweight)
    vi_res = vi - v200_nao*proji
    vp_res = vp - v200_nao*projp

    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:275, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                           my_cmap)

    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)

    caxismin, caxismax, ctrs = -4, 4, np.linspace(-4, 4, 17)

    fig, axs = plt.subplots(3, 2, facecolor='w',
                            edgecolor='k', linewidth=2)

    meshlon, meshlat = np.meshgrid(long, lat)
    meshlon96, meshlat96 = np.meshgrid(lon_n96, lat_n96)
    # subplot 1
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon96, meshlat96)
    plot = m.contourf(x, y, eof_2plot, ctrs,
                      cmap=newcmap1, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[0, 0].annotate('a', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 0], orientation='vertical',
                     spacing='proportional', aspect=50)
    c.set_label(label='NAO (hPa)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 2
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, v200_nao, ctrs/2,
                      cmap=newcmap, vmin=caxismin/2, vmax=caxismax/2,
                      extend='both')
    axs[0, 1].annotate('b', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='v200 NAO (ms$^{-1}$ s.d.$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 3
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vi, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 0].annotate('c', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 0], orientation='vertical',
                     spacing='proportional', aspect=50)
    c.set_label(label='v200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 4
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vi_res, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 1].annotate('d', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 1], orientation='vertical',
                     spacing='proportional', aspect=50)
    c.set_label(label='v200 residual (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 5
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[2, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vp, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[2, 0].annotate('e', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[2, 0], orientation='vertical',
                     spacing='proportional', aspect=50)
    c.set_label(label='v200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 6
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[2, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vp_res, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[2, 1].annotate('f', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[2, 1], orientation='vertical',
                     spacing='proportional', aspect=50)
    c.set_label(label='v200 residual (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)

    for a in range(3):
        for b in range(2):
            m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                        llcrnrlon=0, urcrnrlon=360, resolution='c',
                        ax=axs[a, b])
            ax1 = axs[a, b]
            ax1.set_ylim(-90, 90)
            m.drawparallels(np.arange(-60., 90., 30.),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0., 420., 60.),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-90., 120., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if a == 1:
                polyi = mpl.patches.Polygon([m(45, -20), m(90, -20),
                                             m(90, 00), m(45, 00)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyi)
                #if b == 0:
                   # ax1.set_ylabel('Indian', fontsize=20, labelpad=40)
            if a == 2:
                polyp = mpl.patches.Polygon([m(150, -20), m(195, -20),
                                             m(195, 00), m(150, 00)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyp)

    plt.subplots_adjust(hspace=0.15, wspace=0, top=.99, bottom=0.02,
                        left=.05, right=.95)
    plt.show()
'''


def mapplot90_windres(eof_had, v200_nao, v_sens):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    long = np.concatenate((lon, [360.9375]))
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 360, 193)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    eof_2plot = np.concatenate((eof_had[0, :, 96:], eof_had[0, :, :96]), axis=1)
    eof_2plot = np.concatenate((eof_2plot, np.expand_dims(eof_2plot[:, 0],
                                                          axis=1)), axis=1)
    v200_nao = np.concatenate((v200_nao[0], np.expand_dims(v200_nao[0, :, 0],
                               axis=1)), axis=1)
    vi = np.concatenate((v_sens[0], np.expand_dims(v_sens[0, :, 0], axis=1)), axis=1)
    vp = np.concatenate((v_sens[1], np.expand_dims(v_sens[1, :, 0], axis=1)), axis=1)

    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:275, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                           my_cmap)

    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)

    caxismin, caxismax, ctrs = -4, 4, np.linspace(-4, 4, 17)

    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)

    meshlon, meshlat = np.meshgrid(long, lat)
    meshlon96, meshlat96 = np.meshgrid(lon_n96, lat_n96)
    # subplot 1
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon96, meshlat96)
    plot = m.contourf(x, y, eof_2plot, ctrs,
                      cmap=newcmap1, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[0, 0].annotate('a', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='NAO (hPa)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 2
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, v200_nao, ctrs/2,
                      cmap=newcmap, vmin=caxismin/2, vmax=caxismax/2,
                      extend='both')
    axs[0, 1].annotate('b', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='v200 NAO (ms$^{-1}$ s.d.$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 3
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vi, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 0].annotate('c', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='v200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 4
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vp, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 1].annotate('d', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='v200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)


    for a in range(2):
        for b in range(2):
            m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                        llcrnrlon=0, urcrnrlon=360, resolution='c',
                        ax=axs[a, b])
            ax1 = axs[a, b]
            ax1.set_ylim(-90, 90)
            m.drawparallels(np.arange(-60., 90., 30.),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0., 420., 60.),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-90., 120., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if a == 1 and b == 0:
                polyi = mpl.patches.Polygon([m(45, -20), m(90, -20),
                                             m(90, 00), m(45, 00)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyi)
                #if b == 0:
                   # ax1.set_ylabel('Indian', fontsize=20, labelpad=40)
            if a == 1 and b == 1:
                polyp = mpl.patches.Polygon([m(150, -20), m(195, -20),
                                             m(195, 00), m(150, 00)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyp)

    plt.subplots_adjust(hspace=0.15, wspace=0, top=.99, bottom=0.03,
                        left=.03, right=1)
    plt.show()


def mapplot90_windmslp(mslp, uwnd, uclim200):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    long = np.concatenate((lon, [360.9375]))
    lat1 = np.linspace(90, -90, 145)
    lon1 = np.arange(0, 360, 1.875)
    long1 = np.concatenate((lon1, [360]))

    mslp_i = np.concatenate((mslp[0], np.expand_dims(mslp[0, :, 0],
                                                     axis=1)), axis=1)
    mslp_p = np.concatenate((mslp[1], np.expand_dims(mslp[1, :, 0],
                             axis=1)), axis=1)
    ui = np.concatenate((uwnd[0], np.expand_dims(uwnd[0, :, 0], axis=1)), axis=1)
    up = np.concatenate((uwnd[1], np.expand_dims(uwnd[1, :, 0], axis=1)), axis=1)
    uc = np.concatenate((uclim200[0], np.expand_dims(uclim200[1, :, 0], axis=1)), axis=1)

    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:275, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                           my_cmap)

    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)

    caxismin, caxismax, ctrs = -12, 12, np.linspace(-12, 12, 17)

    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)

    meshlon, meshlat = np.meshgrid(long, lat)
    meshlon1, meshlat1 = np.meshgrid(long1, lat1)
    # subplot 1
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon1, meshlat1)
    plot = m.contourf(x, y, mslp_i/100, ctrs/2,
                      cmap=newcmap1, vmin=caxismin/2, vmax=caxismax/2,
                      extend='both')
    axs[0, 0].annotate('a', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='SLP sensitivity (hPa K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 2
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon1, meshlat1)
    plot = m.contourf(x, y, mslp_p/100, ctrs/2,
                      cmap=newcmap1, vmin=caxismin/2, vmax=caxismax/2,
                      extend='both')
    axs[0, 1].annotate('b', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='SLP sensitivity (hPa K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 3
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, ui, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    m.contour(x, y, uc, [-60, -40, -20, 20, 40, 60], colors='k')
    axs[1, 0].annotate('c', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='u200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 4
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, up, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    m.contour(x, y, uc, [-60, -40, -20, 20, 40, 60], colors='k')
    axs[1, 1].annotate('d', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='u200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)

    for a in range(2):
        for b in range(2):
            m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                        llcrnrlon=0, urcrnrlon=360, resolution='c',
                        ax=axs[a, b])
            ax1 = axs[a, b]
            ax1.set_ylim(-90, 90)
            m.drawparallels(np.arange(-60., 90., 30.),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0., 420., 60.),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-90., 120., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if b == 0:
                polyi = mpl.patches.Polygon([m(45, -20), m(90, -20),
                                             m(90, 00), m(45, 00)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyi)
                #if b == 0:
                   # ax1.set_ylabel('Indian', fontsize=20, labelpad=40)
            if b == 1:
                polyp = mpl.patches.Polygon([m(150, -20), m(195, -20),
                                             m(195, 00), m(150, 00)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyp)

    plt.subplots_adjust(hspace=0.15, wspace=0, top=.99, bottom=0.03,
                        left=.03, right=1)
    plt.show


def mapplot90_windmslp_supp(mslp, uwnd, vwnd, uclim200):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    long = np.concatenate((lon, [360.9375]))
    lat1 = np.linspace(90, -90, 145)
    lon1 = np.arange(0, 360, 1.875)
    long1 = np.concatenate((lon1, [360]))

    mslp_i = np.concatenate((mslp[0], np.expand_dims(mslp[0, :, 0],
                                                     axis=1)), axis=1)
    vi = np.concatenate((vwnd[0], np.expand_dims(vwnd[0, :, 0], axis=1)), axis=1)
    ui = np.concatenate((uwnd[0], np.expand_dims(uwnd[0, :, 0], axis=1)), axis=1)
    uc = np.concatenate((uclim200[0], np.expand_dims(uclim200[1, :, 0], axis=1)), axis=1)

    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:275, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                           my_cmap)

    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)

    caxismin, caxismax, ctrs = -12, 12, np.linspace(-12, 12, 17)

    fig, axs = plt.subplots(3, 1, facecolor='w',
                            edgecolor='k', linewidth=2)

    meshlon, meshlat = np.meshgrid(long, lat)
    meshlon1, meshlat1 = np.meshgrid(long1, lat1)
    # subplot 1
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon1, meshlat1)
    plot = m.contourf(x, y, mslp_i/100, ctrs/2,
                      cmap=newcmap1, vmin=caxismin/2, vmax=caxismax/2,
                      extend='both')
    axs[0].annotate('a', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='SLP sensitivity (hPa K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 2
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[2])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vi, ctrs/3,
                      cmap=newcmap, vmin=caxismin/3, vmax=caxismax/3,
                      extend='both')
    axs[2].annotate('c', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[2], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='v200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 3
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, ui, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    m.contour(x, y, uc, [-60, -40, -20, 20, 40, 60], colors='k')
    axs[1].annotate('b', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='u200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)

    for a in range(3):
        m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=0, urcrnrlon=360, resolution='c',
                    ax=axs[a])
        ax1 = axs[a]
        ax1.set_ylim(-90, 90)
        m.drawparallels(np.arange(-60., 90., 30.),
                        labels=[True, False, False, True], linewidth=0)
        m.drawmeridians(np.arange(0., 420., 60.),
                        labels=[True, False, False, True], linewidth=0)
        ax1.set_yticks(np.arange(-90., 120., 30.))
        ax1.set_xticks(np.arange(0., 390., 30.))
        ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                        direction='out', length=5, width=2)
        ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
        ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
        ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                        direction='out', length=4, width=1)

        polyi = mpl.patches.Polygon([m(300, -5), m(358.5, -5),
                                     m(358.5, 15), m(300, 15)],
                                    linewidth=2, edgecolor='fuchsia',
                                    fill=None)
        ax1.add_patch(polyi)

    plt.subplots_adjust(hspace=0.15, wspace=0, top=.99, bottom=0.03,
                        left=.03, right=1)
    plt.show


def mapplot90_windres_summer(eof_had, v200_nao, v_sens, p_sens):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    long = np.concatenate((lon, [360.9375]))
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 360, 193)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    eof_2plot = np.concatenate((eof_had[1, :, 96:], eof_had[1, :, :96]), axis=1)
    eof_2plot = np.concatenate((eof_2plot, np.expand_dims(eof_2plot[:, 0],
                                                          axis=1)), axis=1)
    v200_nao = np.concatenate((v200_nao[1], np.expand_dims(v200_nao[1, :, 0],
                               axis=1)), axis=1)
    vi = np.concatenate((v_sens[1], np.expand_dims(v_sens[0, :, 0], axis=1)), axis=1)
    vp = np.concatenate((p_sens[1], np.expand_dims(p_sens[0, :, 0], axis=1)), axis=1)

    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:275, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                           my_cmap)

    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)

    caxismin, caxismax, ctrs = -2, 2, np.linspace(-2, 2, 17)

    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)

    meshlon, meshlat = np.meshgrid(long, lat)
    meshlon96, meshlat96 = np.meshgrid(lon_n96, lat_n96)
    # subplot 1
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon96, meshlat96)
    plot = m.contourf(x, y, eof_2plot, ctrs,
                      cmap=newcmap1, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[0, 0].annotate('a', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='NAO (hPa)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 2
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, v200_nao, ctrs/2,
                      cmap=newcmap, vmin=caxismin/2, vmax=caxismax/2,
                      extend='both')
    axs[0, 1].annotate('b', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='v200 NAO (ms$^{-1}$ s.d.$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 3
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon96, meshlat96)
    plot = m.contourf(x, y, vp/100, ctrs,
                      cmap=newcmap1, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 0].annotate('c', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='SLP sensitivity (hPa K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 4
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, vi, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 1].annotate('d', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='v200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)


    for a in range(2):
        for b in range(2):
            m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                        llcrnrlon=0, urcrnrlon=360, resolution='c',
                        ax=axs[a, b])
            ax1 = axs[a, b]
            ax1.set_ylim(-90, 90)
            m.drawparallels(np.arange(-60., 90., 30.),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0., 420., 60.),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-90., 120., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if a == 1 and b == 0:
                polyi = mpl.patches.Polygon([m(100, 10), m(150, 10),
                                             m(150, 30), m(100, 30)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyi)
                #if b == 0:
                   # ax1.set_ylabel('Indian', fontsize=20, labelpad=40)
            if a == 1 and b == 1:
                polyp = mpl.patches.Polygon([m(100, 10), m(150, 10),
                                             m(150, 30), m(100, 30)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyp)

    plt.subplots_adjust(hspace=0.15, wspace=0, top=.99, bottom=0.03,
                        left=.03, right=1)
    plt.show()


def mapplot90_windpr_summer(pr, uwnd):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    long = np.concatenate((lon, [360.9375]))
    lat1 = np.linspace(90, -90, 145)
    lon1 = np.arange(0, 360, 1.875)
    long1 = np.concatenate((lon1, [360]))

    mslp_i = np.concatenate((pr[0], np.expand_dims(pr[0, :, 0],
                                                     axis=1)), axis=1)
    mslp_p = np.concatenate((pr[1], np.expand_dims(pr[1, :, 0],
                             axis=1)), axis=1)
    ui = np.concatenate((uwnd[0], np.expand_dims(uwnd[0, :, 0], axis=1)), axis=1)
    up = np.concatenate((uwnd[1], np.expand_dims(uwnd[1, :, 0], axis=1)), axis=1)

    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:275, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                           my_cmap)

    mycmap = plt.cm.RdBu(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)

    caxismin, caxismax, ctrs = -4, 4, np.linspace(-4, 4, 17)

    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)

    meshlon, meshlat = np.meshgrid(long, lat)
    meshlon1, meshlat1 = np.meshgrid(long1, lat1)
    # subplot 1
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon1, meshlat1)
    plot = m.contourf(x, y, mslp_i, ctrs,
                      cmap=newcmap1, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[0, 0].annotate('a', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='Precip sensitivity (mm day$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 2
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[0, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon1, meshlat1)
    plot = m.contourf(x, y, mslp_p, ctrs,
                      cmap=newcmap1, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[0, 1].annotate('b', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[0, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='Precip sensitivity (mm day$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 3
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 0])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, ui, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 0].annotate('c', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 0], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='u200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)
    # subplot 4
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=axs[1, 1])
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, up, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    axs[1, 1].annotate('d', xy=(-.04, .95), fontname='Arial', fontsize=20,
                       fontweight='bold', xycoords='axes fraction')
    c = fig.colorbar(plot, ax=axs[1, 1], orientation='vertical',
                     spacing='proportional', aspect=50, format='%.0f')
    c.set_label(label='u200 sensitivity (ms$^{-1}$ K$^{-1}$)', size=16)
    for label in c.ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    c.ax.tick_params(labelsize=16)

    for a in range(2):
        for b in range(2):
            m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                        llcrnrlon=0, urcrnrlon=360, resolution='c',
                        ax=axs[a, b])
            ax1 = axs[a, b]
            ax1.set_ylim(-90, 90)
            m.drawparallels(np.arange(-60., 90., 30.),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0., 420., 60.),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-90., 120., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-90., 100., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if b == 0:
                polyi = mpl.patches.Polygon([m(100, 10), m(150, 10),
                                             m(150, 30), m(100, 30)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyi)
                #if b == 0:
                   # ax1.set_ylabel('Indian', fontsize=20, labelpad=40)
            if b == 1:
                polyp = mpl.patches.Polygon([m(120, 10), m(120, 30),
                                             m(205, 30), m(205, 10),
                                             m(120, 10), m(120, -15),
                                             m(205, -15), m(205, 10)],
                                            linewidth=2, edgecolor='fuchsia',
                                            fill=None)
                ax1.add_patch(polyp)

    plt.subplots_adjust(hspace=0.15, wspace=0, top=.99, bottom=0.03,
                        left=.03, right=1)
    plt.show


def scatter_sub(jet_ncep, jet_ncar, jet_era, jet_had):
    def stize(indices):
        # ind_stand = ((indices.transpose((2, 0, 1)) - indices.mean(axis=2))/np.std(indices, axis=2)).transpose(1, 2, 0)
        ind_stand = (indices.transpose((2, 0, 1)) - indices.mean(axis=2)).transpose(1, 2, 0)
        return ind_stand
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1871, 2018)
    jncep = stize(jet_ncep)
    jncar = stize(jet_ncar)
    jera = stize(jet_era)
    #jhad = stize(jet_had)
    # jhad = ((jet_had.transpose((2, 0, 1)) - jet_ncep.mean(axis=2))/np.std(jet_ncep, axis=2)).transpose(1, 2, 0)
    jhad = ((jet_had.transpose((2, 0, 1)) - jet_had.mean(axis=2))).transpose(1, 2, 0)
    for s in range(2):
        for j in range(2):
            for axis in ['top', 'bottom', 'left', 'right']:
                axs[s, 1-j].spines[axis].set_linewidth(2)
            cc = np.corrcoef(jncep[s, j], jhad[s, j, 109:])[0, 1]
            cc1 = np.corrcoef(jncar[s, j], jhad[s, j, 1:142])[0, 1]
            cc2 = np.corrcoef(jera[s, j], jhad[s, j, 30:140])[0, 1]
            axs[s, 1-j].plot(time[109:], jncep[s, j], color='#0247FE', label='NCEP2: r=%.2f' % cc)
            axs[s, 1-j].plot(time[1:142], jncar[s, j], color='#700CBC', label='NOAA-20C: r=%.2f' % cc1)
            axs[s, 1-j].plot(time[30:140], jera[s, j], color='#66B032', label='ERA-20C: r=%.2f' % cc2)
            axs[s, 1-j].plot(time, jhad[s, j], color='r', label='Reconstructed')
            axs[s, 1-j].set_xlim([1871, 2017])
            axs[s, 1-j].set_xlabel('Year', fontsize=16)
            axs[s, 1-j].legend(loc=2, fontsize=16)
            if j == 0:
                axs[s, 1-j].set_ylabel('Jet speed (ms$\mathregular{^{-1}}$)', fontsize=16)
                axs[s, 1-j].set_ylim([-2.5, 2.5])
                if s == 0:
                    axs[s, 1-j].set_title('Speed', fontsize=16)
                    axs[s, 1-j].set_ylim([-4, 4])
            if j == 1:
                axs[s, 1-j].set_ylabel('Jet position (deg)', fontsize=16, labelpad=10)
                axs[s, 1-j].set_ylim([-6.5, 6.5])
                if s == 0:
                    axs[s, 1-j].set_ylabel('Jet position (deg)', fontsize=16, labelpad=-10)
                    axs[s, 1-j].set_title('Latitude', fontsize=16)
                    axs[s, 1-j].set_ylim([-10.5, 10.5])
            axs[s, 1-j].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-65, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)

def scatter_sub5(jet_ncep, jet_ncar, jet_era, jet_had):
    def stize(indices):
        ind_stand = ((indices.transpose((2, 0, 1)) - indices.mean(axis=2))/np.std(indices, axis=2)).transpose(1, 2, 0)
        return ind_stand
    def lanczos(field, wi=51, co=0.0125, dim=1):
        def lweights(window, cutoff):
            """
            Calculate weights for a low pass Lanczos filter
            and then applies to time series
    
            Parameters
            ----------
    
            window: int
                The length of the filter window.
            cutoff: float
                The cutoff frequency in inverse time steps.
    
            Returns
            -------
            w: array
                filter weights
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
            w = w / np.sum(w)
            return w

        l_filter = lweights(wi, co)

        if dim == 1:
            filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1))
            for i in range(len(field)-len(l_filter)+1):
                filtered[i] = np.sum(field[i:i+len(l_filter)] * l_filter)
        if dim == 2:
            filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1,
                                np.ma.size(field, axis=1)))
            for i in range(len(field)-len(l_filter)+1):
                filtered[i] = np.sum((field[i:i+len(l_filter)].transpose() *
                                      l_filter).transpose(), axis=0)
        return filtered
    jet_had5 = np.zeros((2, 2, 137))
    jet_ncep5 = np.zeros((2, 2, 28))
    jet_ncar5 = np.zeros((2, 2, 131))
    jet_era5 = np.zeros((2, 2, 100))
    for i in range(2):
        for j in range(2):
            jet_had5[i, j] = lanczos(jet_had[i, j], 11, 1/5)
            jet_ncep5[i, j] = lanczos(jet_ncep[i, j], 11, 1/5)
            jet_ncar5[i, j] = lanczos(jet_ncar[i, j], 11, 1/5)
            jet_era5[i, j] = lanczos(jet_era[i, j], 11, 1/5)
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1876, 2013)
    jncep = stize(jet_ncep5)
    jncar = stize(jet_ncar5)
    jera = stize(jet_era5)
    jhad = stize(jet_had5)
    for s in range(2):
        for j in range(2):
            for axis in ['top', 'bottom', 'left', 'right']:
                axs[s, 1-j].spines[axis].set_linewidth(2)
            cc = np.corrcoef(jncep[s, j], jhad[s, j, 109:])[0, 1]
            cc1 = np.corrcoef(jncar[s, j], jhad[s, j, 1:132])[0, 1]
            cc2 = np.corrcoef(jera[s, j], jhad[s, j, 30:130])[0, 1]
            axs[s, 1-j].plot(time[109:], jncep[s, j], color='#0247FE', label='NCEP2: r=%.2f' % cc)
            axs[s, 1-j].plot(time[1:132], jncar[s, j], color='#700CBC', label='NOAA-20C: r=%.2f' % cc1)
            axs[s, 1-j].plot(time[30:130], jera[s, j], color='#66B032', label='ERA-20C: r=%.2f' % cc2)
            axs[s, 1-j].plot(time, jhad[s, j], color='r', label='Reconstructed')
            axs[s, 1-j].set_ylim([-3, 3])
            axs[s, 1-j].set_xlim([1876, 2013])
            axs[s, 1-j].set_xlabel('Year', fontsize=16)
            axs[s, 1-j].legend(loc=2, fontsize=16)
            if j == 0:
                axs[s, 1-j].set_ylabel('Jet speed (ms$\mathregular{^{-1}}$)', fontsize=16)
                if s == 0:
                    axs[s, 1-j].set_title('Speed', fontsize=16)
            if j == 1:
                axs[s, 1-j].set_ylabel('Jet position (deg)', fontsize=16)
                if s == 0:
                    axs[s, 1-j].set_title('Latitude', fontsize=16)
            axs[s, 1-j].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def moving_plot30jet(jet_ncep, jet_ncar, jet_era, jet_had):
    def moving_corr30(jet_ncep, jet_ncar, jet_era, jet_rec):
        ncep = np.zeros((2, 2, 8))
        ncar = np.zeros((2, 2, 111))
        era = np.zeros((2, 2, 80))
        for s in range(2):
            for t in range(2):
                jet_rec_s = (jet_rec[s, t] - np.mean(jet_rec[s, t])) / np.std(jet_rec[s, t])
                for i in range(8):
                    ncep[s, t, i] = np.corrcoef(jet_ncep[s, t, i:i+31], jet_rec_s[109+i:109+i+31])[0, 1]
                for i in range(111):
                    ncar[s, t, i] = np.corrcoef(jet_ncar[s, t, i:i+31], jet_rec_s[1+i:1+i+31])[0, 1]
                for i in range(80):
                    era[s, t, i] = np.corrcoef(jet_era[s, t, i:i+31], jet_rec_s[30+i:30+i+31])[0, 1]
        return ncep, ncar, era

    import matplotlib.pyplot as plt
    m_ncep, m_ncar, m_era = moving_corr30(jet_ncep, jet_ncar, jet_era, jet_had)
    fig, axs = plt.subplots(2, 2, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1886, 2003)
    for s in range(2):
        for j in range(2):
            for axis in ['top', 'bottom', 'left', 'right']:
                axs[s, 1-j].spines[axis].set_linewidth(2)
            axs[s, 1-j].plot(time[109:], m_ncep[s, j], color='#0247FE', label='NCEP2')
            axs[s, 1-j].plot(time[1:112], m_ncar[s, j], color='#700CBC', label='NOAA-20C')
            axs[s, 1-j].plot(time[30:110], m_era[s, j], color='#66B032', label='ERA-20C')
            axs[s, 1-j].axhline(rsig(30), ls='--', color='r')
            axs[s, 1-j].axhline(-rsig(30), ls='--', color='r')
            axs[s, 1-j].axhline(color='k')  
            axs[s, 1-j].set_ylim([-.4, .6])
            axs[s, 1-j].set_xlim([1886, 2003])
            axs[s, 1-j].set_xlabel('Year', fontsize=16)
            axs[s, 1-j].legend(loc=2, fontsize=16)
            if j == 0:
                if s == 0:
                    axs[s, 1-j].set_title('Speed', fontsize=16)
            if j == 1:
                if s == 0:
                    axs[s, 1-j].set_title('Latitude', fontsize=16)
            axs[s, 1-j].tick_params(axis='both', which='major', labelsize=16)
    #plt.suptitle('Jet indices 31 year moving average correlation', fontsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def scatter_subnao(nao_ncep, nao_ncar, nao_era, nao_rec):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1871, 2018)
    for s in range(2):
        for axis in ['top', 'bottom', 'left', 'right']:
            axs[s].spines[axis].set_linewidth(2)
        nao_rec_s = (nao_rec[s]) - np.mean(nao_rec[s])#) / np.std(nao_rec[s])
        #nao_rec_s = lanczos(nao_rec_s, 31, 1/30, hl='low')
        cc = np.corrcoef(nao_ncep[s], nao_rec_s[109:])[0, 1]
        cc1 = np.corrcoef(nao_ncar[s], nao_rec_s[1:142])[0, 1]
        cc2 = np.corrcoef(nao_era[s], nao_rec_s[30:140])[0, 1]
        axs[s].plot(time[109:], nao_ncep[s], color='#0247FE', label='NCEP2: r=%.2f' % cc)
        axs[s].plot(time[1:142], nao_ncar[s], color='#700CBC', label='NOAA-20C: r=%.2f' % cc1)
        axs[s].plot(time[30:140], nao_era[s], color='#66B032', label='ERA-20C: r=%.2f' % cc2)
        axs[s].plot(time, nao_rec_s, color='r', label='Reconstructed')
        axs[s].axhline(color='k')
        axs[s].set_ylim([-3, 3])
        axs[s].set_xlim([1871, 2017])
        axs[s].set_xlabel('Year', fontsize=16)
        axs[s].legend(loc=2, fontsize=16)
        if s == 0:
            axs[s].set_title('NAO', fontsize=16)
        axs[s].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs, ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def scatter_subnao_att(nao_ncar, nao_glo, nao_atl, nao_ind, nao_pac, winter='yes'):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(5, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1871, 2018)
    if winter == 'yes':
        s = 0
        l = 1
        p = 1.5
    else:
        s = 1
        l = .3
        p = .45
    for i in range(5):
        for axis in ['top', 'bottom', 'left', 'right']:
            axs[i].spines[axis].set_linewidth(2)
            if i != 0:
                axs[i].spines['top'].set_visible(False)
            if i != 4:
                axs[i].spines['bottom'].set_visible(False)
    nao_glo1 = (nao_glo[s]) - np.mean(nao_glo[s])
    nao_pac1 = (nao_pac[s]) - np.mean(nao_pac[s])
    nao_ind1 = (nao_ind[s]) - np.mean(nao_ind[s])
    nao_atl1 = (nao_atl[s]) - np.mean(nao_atl[s])
    cc = np.corrcoef(nao_ncar[s], nao_glo1[1:-5])[0, 1]
    cc1 = np.corrcoef(nao_ncar[s], nao_ind1[1:-5])[0, 1]
    cc2 = np.corrcoef(nao_ncar[s], nao_pac1[1:-5])[0, 1]
    cc3 = np.corrcoef(nao_ncar[s], nao_atl1[1:-5])[0, 1]

    axs[0].plot(time[1:142], nao_ncar[s], color='r')
    if winter == 'no':
        line = np.zeros((70))
        line[:] = 2.4
        axs[0].plot(time[54:124], line, lw=2, color='blue')
        axs[0].plot([1925, 1925], [2.2, 2.6], lw=2, color='blue')
        axs[0].plot([1994, 1994], [2.2, 2.6], lw=2, color='blue')
    axs[1].plot(time, nao_glo1, color='r', label='Global: r=%.2f' % cc)
    axs[2].plot(time, nao_atl1, color='r', label='Atl: r=%.2f' % cc3)
    axs[3].plot(time, nao_ind1, color='r', label='Ind: r=%.2f' % cc1)
    axs[2].plot(time, nao_ind1-3, color='r')
    axs[4].plot(time, nao_pac1, color='r', label='Pac: r=%.2f' % cc2)
    if winter != 'yes':
        axs[0].plot(time, nao_glo1*20/3-6, color='r')
    axs[0].set_ylim([-3, 3])
    axs[1].set_yticks([-l, 0, l])
    axs[2].set_yticks([-l, 0, l])
    axs[3].set_yticks([-l, 0, l])
    axs[4].set_yticks([-l, 0, l])
    axs[0].axhline(color='k')
    axs[1].axhline(color='k')
    axs[2].axhline(color='k')
    axs[3].axhline(color='k')
    axs[4].axhline(color='k')
    axs[1].set_ylim([-p, p])
    axs[2].set_ylim([-p, p])
    axs[3].set_ylim([-p, p])
    axs[4].set_ylim([-p, p])
    for i in range(1880, 2020, 10):
        for j in range(5):
            axs[j].axvline(i, ls='--', color='k')
    for i in range(5):
        axs[i].set_xlim([1871, 2017])
        axs[i].tick_params(axis='both', which='major', labelsize=16)
    axs[4].set_xlabel('Year', fontsize=16)
    #axs[0].legend(loc=2, fontsize=16)
    rows = ['NOAA-20thC\n', 'Global\n(0.38, 0.38)', 'Atl\n(0.27, 0.31)', 'Ind\n(0.01, 0.58)', 'Pac\n(0.16, 0.32)']
    #rows = ['NOAA-20thC\n', 'Global\n(0.37, 0.16)', 'Atl\n(0.10, 0.06)', 'Ind\n(-0.01, 0.07)', 'Pac\n(0.35, 0.15)']
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-65, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='center', rotation=90, fontsize=20)

    plt.subplots_adjust(hspace=0, wspace=0.15, top=.99, bottom=0.075, left=.075,
                        right=.99)


def pacific_std(m_ncar, std_pac_warm):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 1, figsize=(3.8, 2.6),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1870+15, 2018-15)
    for axis in ['top', 'bottom', 'left', 'right']:
        axs.spines[axis].set_linewidth(2)

    m_ncar1 = (m_ncar[1] - np.mean(m_ncar[1]))/np.std(m_ncar[1])
    std1 = (std_pac_warm - np.mean(std_pac_warm))/np.std(std_pac_warm)
    cc = np.corrcoef(m_ncar1, std_pac_warm[2:-5])[0, 1]

    axs.plot(time[2:-5], m_ncar1, color='k', label='Reconstruction skill')
    axs.plot(time, std1, color='r', label='Pacific SST s.d.')
    axs.plot(time, std1-100, color='white', label='Correlation: %.2f' % cc)
    axs.plot(time[2:-5], yax, color='blue')
    axs.set_yticks([-2, -1, 0, 1, 2])
    axs.set_ylim([-2, 2.5])
    axs.set_xlim([1870+15+2, 2017-15-5])
    axs.tick_params(axis='both', which='major', labelsize=11)
    axs.set_xlabel('Year', fontsize=8)
    #plt.text(1890, 2,'Correlation: %.2f' % cc, fontsize=16)
    axs.legend(loc=4, fontsize=8)
    plt.subplots_adjust(hspace=0, wspace=0.15, top=.975, bottom=0.1,
                        left=.1, right=.975)


def scatter_subnao_lanczos(nao_ncep, nao_ncar, nao_era, nao_rec, hil='low'):
    import matplotlib.pyplot as plt
    from rpm_2 import lanczos
    fig, axs = plt.subplots(2, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1871+15, 2018-15)
    for s in range(2):
        for axis in ['top', 'bottom', 'left', 'right']:
            axs[s].spines[axis].set_linewidth(2)
        nao_rec_s = (nao_rec[s])# - np.mean(nao_rec[s])) / np.std(nao_rec[s])
        nao_rec_s = lanczos(nao_rec_s, 31, 1/10, hl=hil)
        nao_ncep1 = lanczos(nao_ncep[s], 31, 1/10, hl=hil)
        nao_ncar1 = lanczos(nao_ncar[s], 31, 1/10, hl=hil)
        nao_era1 = lanczos(nao_era[s], 31, 1/10, hl=hil)
        cc = np.corrcoef(nao_ncep1, nao_rec_s[109:])[0, 1]
        cc1 = np.corrcoef(nao_ncar1, nao_rec_s[1:142-30])[0, 1]
        cc2 = np.corrcoef(nao_era1, nao_rec_s[30:140-30])[0, 1]
        axs[s].plot(time[109:], nao_ncep1, color='#0247FE', label='NCEP2: r=%.2f' % cc)
        axs[s].plot(time[1:142-30], nao_ncar1, color='#700CBC', label='NOAA-20C: r=%.2f' % cc1)
        axs[s].plot(time[30:140-30], nao_era1, color='#66B032', label='ERA-20C: r=%.2f' % cc2)
        axs[s].plot(time, nao_rec_s, color='r', label='Reconstructed')
        if hil == 'low':
            axs[s].set_ylim([-2, 2])
        else:
            axs[s].set_ylim([-3, 3])
        axs[s].set_xlim([1871+15, 2017-15])
        axs[s].set_xlabel('Year', fontsize=16)
        axs[s].legend(loc=2, fontsize=16)
        if s == 0:
            axs[s].set_title('NAO', fontsize=16)
        axs[s].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs, ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def scatter_subnao5(nao_ncep, nao_ncar, nao_era, nao_had):
    def lanczos(field, wi=51, co=0.0125, dim=1, hl='low'):
        def lweights(window, cutoff, hilo='low'):
            """
            Calculate weights for a low pass Lanczos filter
            and then applies to time series
    
            Parameters
            ----------
    
            window: int
                The length of the filter window.
            cutoff: float
                The cutoff frequency in inverse time steps.
    
            Returns
            -------
            w: array
                filter weights
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
            w = w / np.sum(w)
            if hilo == 'high':
                w = w*-1
                w[order-1] += 1
            return w
    
        l_filter = lweights(wi, co, hilo=hl)
    
        if dim == 1:
            filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1))
            for i in range(len(field)-len(l_filter)+1):
                filtered[i] = np.sum(field[i:i+len(l_filter)] * l_filter)
        if dim == 2:
            filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1,
                                np.ma.size(field, axis=1)))
            for i in range(len(field)-len(l_filter)+1):
                filtered[i] = np.sum((field[i:i+len(l_filter)].transpose() *
                                      l_filter).transpose(), axis=0)
        if dim == 3:
            filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1,
                                np.ma.size(field, axis=1),
                                np.ma.size(field, axis=2)))
            for i in range(len(field)-len(l_filter)+1):
                filtered[i] = np.sum((field[i:i+len(l_filter)].transpose((1, 2, 0)) *
                                      l_filter).transpose((2, 0, 1)), axis=0)
        return filtered
    nao_had5 = np.zeros((2, 117))
    nao_ncep5 = np.zeros((2, 8))
    nao_ncar5 = np.zeros((2, 111))
    nao_era5 = np.zeros((2, 80))
    for i in range(2):
        nao_had5[i] = lanczos(nao_had[i], 31, 1/30, hl='high')
        nao_ncep5[i] = lanczos(nao_ncep[i], 31, 1/30, hl='high')
        nao_ncar5[i] = lanczos(nao_ncar[i], 31, 1/30, hl='high')
        nao_era5[i] = lanczos(nao_era[i], 31, 1/30, hl='high')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1886, 2003)
    for s in range(2):
        for axis in ['top', 'bottom', 'left', 'right']:
            axs[s].spines[axis].set_linewidth(2)
        #nao_rec_s = (nao_rec[s] - np.mean(nao_rec[s])) / np.std(nao_rec[s])
        cc = np.corrcoef(nao_ncep5[s], nao_had5[s, 109:109+8])[0, 1]
        cc1 = np.corrcoef(nao_ncar5[s], nao_had5[s, 1:-5])[0, 1]
        cc2 = np.corrcoef(nao_era5[s], nao_had5[s, 30:-7])[0, 1]
        axs[s].plot(time[109:109+8], nao_ncep5[s], color='#0247FE', label='NCEP2: r=%.2f' % cc)
        axs[s].plot(time[1:-5], nao_ncar5[s], color='#700CBC', label='NOAA-20C: r=%.2f' % cc1)
        axs[s].plot(time[30:-7], nao_era5[s], color='#66B032', label='ERA-20C: r=%.2f' % cc2)
        axs[s].plot(time, nao_had5[s], color='r', label='Reconstructed')
        axs[s].set_ylim([-3, 3])
        axs[s].set_xlim([1886, 2002])
        axs[s].set_xlabel('Year', fontsize=16)
        axs[s].legend(loc=2, fontsize=16)
        if s == 0:
            axs[s].set_title('NAO', fontsize=16)
        axs[s].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs, ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def moving_plot30(nao_ncep, nao_ncar, nao_era, nao_had):
    def moving_corr30(nao_ncep, nao_ncar, nao_era, nao_rec):
        ncep = np.zeros((2, 8))
        ncar = np.zeros((2, 111))
        #ncar_neff = np.zeros((2, 111))
        era = np.zeros((2, 80))
        for s in range(2):
            nao_rec_s = (nao_rec[s] - np.mean(nao_rec[s])) / np.std(nao_rec[s])
            for i in range(8):
                ncep[s, i] = np.corrcoef(nao_ncep[s, i:i+31], nao_rec_s[109+i:109+i+31])[0, 1]
            for i in range(111):
                ncar[s, i] = np.corrcoef(nao_ncar[s, i:i+31], nao_rec_s[1+i:1+i+31])[0, 1]
                #rho1 = np.corrcoef(nao_ncar[s, i+1:i+31], nao_ncar[s, i:i+30])[0, 1]
                #rho2 = np.corrcoef(nao_rec_s[2+i:1+i+31], nao_rec_s[1+i:1+i+30])[0, 1]
                #ncar_neff[s, i] = rsig(30 * (1-rho1*rho2)/(1+rho1*rho2))
                #ncar[s, i] = np.sum((nao_rec_s[1+i:1+i+31])*(nao_ncar[s, i:i+31]))/(np.sqrt(np.sum((nao_rec_s[1+i:1+i+31])**2))*np.sqrt(np.sum((nao_ncar[s, i:i+31])**2)))
            for i in range(80):
                era[s, i] = np.corrcoef(nao_era[s, i:i+31], nao_rec_s[30+i:30+i+31])[0, 1]
        return ncep, ncar, era#, ncar_neff

    m_ncep, m_ncar, m_era = moving_corr30(nao_ncep, nao_ncar, nao_era, nao_had)

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1886, 2003)
    for s in range(2):
        for axis in ['top', 'bottom', 'left', 'right']:
            axs[s].spines[axis].set_linewidth(2)
        axs[s].plot(time[109:], m_ncep[s], color='#0247FE', label='NCEP2')
        axs[s].plot(time[1:112], m_ncar[s], color='#700CBC', label='NOAA-20C')
        #axs[s].plot(time[1:112], ncar_neff[s], color='r', label='NOAA-20C')
        axs[s].plot(time[30:110], m_era[s], color='#66B032', label='ERA-20C')
        axs[s].axhline(rsig(29.407583564578303), ls='--', color='r')
        axs[s].axhline(-rsig(29.407583564578303), ls='--', color='r')
        axs[s].axhline(color='k')  
        axs[s].set_ylim([-.4, .7])
        axs[s].set_xlim([1886, 2003])
        axs[s].set_xlabel('Year', fontsize=16)
        axs[s].legend(loc=2, fontsize=16)
        axs[s].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs, ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def moving_plot30_5(nao_ncar, nao_era, nao_had):
    def moving_corr30_5(nao_ncar, nao_era, nao_rec):
        ncar = np.zeros((2, 81))
        era = np.zeros((2, 50))
        for s in range(2):
            nao_rec_s = (nao_rec[s] - np.mean(nao_rec[s])) / np.std(nao_rec[s])
            for i in range(81):
                ncar[s, i] = np.corrcoef(nao_ncar[s, i:i+31], nao_rec_s[1+i:1+i+31])[0, 1]
            for i in range(50):
                era[s, i] = np.corrcoef(nao_era[s, i:i+31], nao_rec_s[30+i:30+i+31])[0, 1]
        return ncar, era

    m_ncar, m_era = moving_corr30_5(nao_ncar, nao_era, nao_had)

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1901, 1988)
    for s in range(2):
        axs[s].plot(time[1:-5], m_ncar[s], color='#700CBC', label='NOAA-20C')
        axs[s].plot(time[30:-7], m_era[s], color='#66B032', label='ERA-20C')
        axs[s].set_ylim([-1, 1])
        axs[s].set_xlim([1901, 1988])
        axs[s].set_xlabel('Year', fontsize=16)
        axs[s].legend(loc=2, fontsize=16)
        if s == 0:
            axs[s].set_title('NAO 31 year moving average correlation', fontsize=16)
        axs[s].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs, ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def scatter_ext(jet_ncep, jet_ncar, pr_djf, jet_had):
    def stize(indices):
        ind_stand = (indices - indices.mean())
        return ind_stand
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1871, 2018)
    jncep = stize(jet_ncep)
    jncar = stize(jet_ncar)
    jera = stize(pr_djf)
    #jhad = stize(jet_had)
    jhad = (jet_had - jet_had.mean())
    for axis in ['top', 'bottom', 'left', 'right']:
        axs.spines[axis].set_linewidth(2)
    cc = np.corrcoef(jncep, jhad[109:])[0, 1]
    cc1 = np.corrcoef(jncar, jhad[1:142])[0, 1]
    cc2 = np.corrcoef(jera, jhad[40:])[0, 1]
    lns1 = axs.plot(time[109:], jncep, color='#0247FE', label='NCEP2: r=%.2f' % cc)
    lns2 = axs.plot(time[1:142], jncar, color='#700CBC', label='NOAA-20C: r=%.2f' % cc1)
    lns4 = axs.plot(time, jhad, color='r', label='Reconstructed')
    axs.set_ylim([-12, 12])
    axs.set_xlabel('Year', fontsize=16)
    axs.set_ylabel('Jet extension anomaly (ms$^{-1}$)', fontsize=16)
    axs.tick_params(axis='both', which='major', labelsize=16)
    ax2 = axs.twinx()
    ax2.set_ylabel('Precipitation (mm season$^{-1}$)', fontsize=16, rotation=-90)
    lns3 = ax2.plot(time[40:], jera ,color='#66B032', label='Precip: r=%.2f' % cc2)
    ax2.tick_params(axis='both', which='major', labelsize=16)
    lns = lns1 + lns2 + lns3 + lns4
    labs = [l.get_label() for l in lns]
    axs.legend(lns, labs, loc=2, fontsize=16)
    axs.set_xlim([1871, 2017])
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.92)


def lagged_plots(c_ncep, c_ncar, c_era, ws='winter'):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    if ws == 'winter':
        xlab = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 
                'Sep', 'Oct', 'Nov', 'Dec']
        #xlab = ['JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS', 'ASO', 
            #    'SON', 'OND', 'NDJ', 'DJF']
        se = 0
    else:
        xlab = ['Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb',
                'Mar', 'Apr', 'May', 'Jun']
        #xlab = ['JAS', 'ASO', 'SON', 'OND', 'NDJ', 'DJF', 'JFM', 'FMA',
            #    'MAM', 'AMJ', 'MJJ', 'JJA']
        se = 1
    for s in range(2):
        for j in range(2):
            for axis in ['top', 'bottom', 'left', 'right']:
                axs[s, j].spines[axis].set_linewidth(2)
            axs[s, j].plot(c_ncep[se, j+2*s], color='#0247FE', label='NCEP2')
            axs[s, j].plot(c_ncar[se, j+2*s], color='#700CBC', label='NOAA-20C')
            axs[s, j].plot(c_era[se, j+2*s], color='#66B032', label='ERA-20C')
            axs[s, j].set_ylim([-.2, .4])
            axs[s, j].set_ylabel('Correlation', fontsize=16)
            #axs[s, j].set_xlim([1876, 2013])
            axs[s, j].legend(loc=2, fontsize=16)
            if s == 0 and j == 0:
                axs[s, j].set_title('Global', fontsize=16)
            if s == 0 and j == 1:
                axs[s, j].set_title('Atlantic', fontsize=16)
            if s == 1 and j == 0:
                axs[s, j].set_title('Indian', fontsize=16)
            if s == 1 and j == 1:
                axs[s, j].set_title('Pacific', fontsize=16)
            axs[s, j].axhline(rsig(35.7749033607328), ls='--', color='#0247FE')
            axs[s, j].axhline(rsig(126.93103340932541), ls='--', color='#700CBC')
            axs[s, j].axhline(rsig(88.547504045570435), ls='--', color='#66B032')
            axs[s, j].axhline(-rsig(35.7749033607328), ls='--', color='#0247FE')
            axs[s, j].axhline(-rsig(126.93103340932541), ls='--', color='#700CBC')
            axs[s, j].axhline(-rsig(88.547504045570435), ls='--', color='#66B032')
            axs[s, j].axhline(color='k')  
            axs[s, j].set_xticks(np.arange(0, 12))
            axs[s, j].set_xticklabels(xlab)
            axs[s, j].tick_params(axis='both', which='major', labelsize=16)

    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def lagged_plots2(c_ncep, c_ncar, c_era, ch_ncep, ch_ncar, ch_era, ws='winter'):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    if ws == 'winter':
        xlab = ['Jun', 'Jul', 'Aug', 
                'Sep', 'Oct', 'Nov', 'Dec', 'DJF']
        #xlab = ['JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS', 'ASO', 
            #    'SON', 'OND', 'NDJ', 'DJF']
        se = 0
    else:
        xlab = ['Dec', 'Jan', 'Feb',
                'Mar', 'Apr', 'May', 'Jun', 'JJA']
        #xlab = ['JAS', 'ASO', 'SON', 'OND', 'NDJ', 'DJF', 'JFM', 'FMA',
            #    'MAM', 'AMJ', 'MJJ', 'JJA']
        se = 1
    for s in range(2):
        for j in range(2):
            for axis in ['top', 'bottom', 'left', 'right']:
                axs[s, j].spines[axis].set_linewidth(2)
            axs[s, j].plot(c_ncep[se, j+2*s, 5:], color='#0247FE', label='NCEP2')
            axs[s, j].plot(c_ncar[se, j+2*s, 5:], color='#700CBC', label='NOAA-20C')
            axs[s, j].plot(c_era[se, j+2*s, 5:], color='#66B032', label='ERA-20C')
            axs[s, j].scatter(7, ch_ncep[se, j+2*s], marker='x', color='#0247FE')
            axs[s, j].scatter(7, ch_ncar[se, j+2*s], marker='x', color='#700CBC')
            axs[s, j].scatter(7, ch_era[se, j+2*s], marker='x', color='#66B032')
            axs[s, j].set_ylim([-.325, .425])
            axs[s, j].set_ylabel('Correlation', fontsize=16)
            axs[s, j].legend(loc=2, fontsize=16)
            if s == 0 and j == 0:
                axs[s, j].set_title('Global', fontsize=16)
            if s == 0 and j == 1:
                axs[s, j].set_title('Atlantic', fontsize=16)
            if s == 1 and j == 0:
                axs[s, j].set_title('Indian', fontsize=16)
            if s == 1 and j == 1:
                axs[s, j].set_title('Pacific', fontsize=16)
            axs[s, j].axhline(rsig(35.7749033607328), ls='--', color='#0247FE')
            axs[s, j].axhline(rsig(126.93103340932541), ls='--', color='#700CBC')
            axs[s, j].axhline(rsig(88.547504045570435), ls='--', color='#66B032')
            axs[s, j].axhline(-rsig(35.7749033607328), ls='--', color='#0247FE')
            axs[s, j].axhline(-rsig(126.93103340932541), ls='--', color='#700CBC')
            axs[s, j].axhline(-rsig(88.547504045570435), ls='--', color='#66B032')
            axs[s, j].axhline(color='k')  
            axs[s, j].set_xticks(np.arange(0, 8))
            axs[s, j].set_xticklabels(xlab)
            axs[s, j].set_xlim([-.25, 7.25])
            axs[s, j].tick_params(axis='both', which='major', labelsize=16)

    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def rsig(n, sig=0.05, twotail='yes'):
    from scipy import stats
    if twotail == 'yes':
        tt = 2
    else:
        tt = 1
    t = stats.norm.ppf(1-sig/tt)
    r = np.sqrt(t**2/(t**2+n-2))
    return r


def mapplot60_je(plotdata_sig, plotdata, mx=2.5, mask='yes', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    long = np.concatenate((lon, [360]))
    fig, ax1 = plt.subplots(1, 1, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    plotdata1 = np.concatenate((plotdata_sig,
                                np.expand_dims(plotdata[:, 0],
                                axis=1)), axis=1)
    plotdata2 = np.concatenate((plotdata,
                                np.expand_dims(plotdata[:, 0],
                                axis=1)), axis=1)
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                             axis=1)
        plotdata1 = np.ma.masked_array(plotdata1, lsm)
        plotdata2 = np.ma.masked_array(plotdata2, lsm)
    meshlon, meshlat = np.meshgrid(long, lat)
    m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
    m.drawcoastlines()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    # my_cmap[238:275, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                           my_cmap)
    caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
    plot = m.contourf(x, y, plotdata1, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    m.contour(x, y, plotdata2, ctrs, colors='k')
    ax1.set_ylim(-60, 60)
    m.drawparallels(np.arange(-60., 90., 30.),
                    labels=[True, False, False, True], linewidth=0)
    m.drawmeridians(np.arange(0., 390., 30.),
                    labels=[True, False, False, True], linewidth=0)
    ax1.set_yticks(np.arange(-60., 90., 30.))
    ax1.set_xticks(np.arange(0., 390., 30.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
    ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    poly = mpl.patches.Polygon([m(long[160], lat[60]), m(359.5, lat[60]),
                                m(359.5, 59.5), m(long[160], 59.5)],
                               linewidth=2, edgecolor='fuchsia', fill=None)
    plt.gca().add_patch(poly)
    poly = mpl.patches.Polygon([m(30, -30), m(105, -30),
                                m(105, 20), m(30, 20)],
                               linewidth=2, edgecolor='fuchsia', fill=None)
    plt.gca().add_patch(poly)
    poly = mpl.patches.Polygon([m(105, -20), m(240, -20),
                                m(240, 30), m(105, 30)],
                               linewidth=2, edgecolor='fuchsia', fill=None)
    plt.gca().add_patch(poly)

    c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                     spacing='proportional', aspect=50)
    c.set_label(label='Extension sensitivity (ms$^{-1} [\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    #c.set_label(label='Circulation index sensitivity (ms$^{-1}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    plt.subplots_adjust(hspace=0, wspace=0.1, top=.95, bottom=0.05, left=.05,
                        right=.95)
    plt.show()


def animate(data, winter='yes', save='no',
            location='/home/bakerh/Downloads/new.mp4'):
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
    import matplotlib.animation as animation
    from mpl_toolkits.basemap import Basemap, addcyclic
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    if winter == 'yes':
        colormax = 5
        colormin = -colormax
        ws = 0
    else:
        colormax = 2.5
        colormin = -colormax
        ws = 1

    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)

    fig = plt.figure(figsize=(16, 10.5), facecolor='w', edgecolor='k', linewidth=2)
    pdata = data[:, ws]
    pdata, lon1 = addcyclic(pdata, lon)
    meshlon, meshlat = np.meshgrid(lon1, lat)
    ctrs = np.linspace(colormin, colormax, 17)

    def updatefig(i):
        fig.clear()
        pdata1 = pdata[i]
        if pdata1[30, 80] < 0: # 44 for winter
            pdata1 *= -1
        ax1 = fig.add_subplot(1, 1, 1)
        m = Basemap(width=10000000, height=7000000,
                    resolution='c', projection='aea',
                    lat_1=40., lat_2=60, lon_0=-25, lat_0=50, ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        plot = m.contourf(x, y, pdata1, ctrs,
                          cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                          extend='both')
        parallels = m.drawparallels(np.arange(-90., 75., 15.))
        meridians = m.drawmeridians(np.arange(-180., 181., 30))
        m.drawparallels(parallels, labels=[True, False, False, False])
        m.drawmeridians(meridians, labels=[False, False, False, True])
        ax1.set_title(str(i+1897), fontname='Arial', fontsize=16)
        cbar_ax = fig.add_axes([0.3, 0.06, 0.4, 0.01])
        b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                         orientation='horizontal', extend='max', format='%.2f')
        b.set_label(label='NAO (hPa per SD)', size=16,
                    fontsize=16, fontname='Arial', labelpad=0)
        cl = plt.getp(cbar_ax, 'xmajorticklabels')
        plt.setp(cl, fontname='Arial', fontsize=16)
        for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
        cbar_ax.tick_params(labelsize=16)
        plt.draw()

    anim = animation.FuncAnimation(fig, updatefig, len(data))
    if save == 'yes':
        anim.save(location, codec='mpeg4', bitrate=8000, dpi=300)
    return anim


def fig_schaller(j, j_na, j_in, j_pa):
    import matplotlib.pyplot as plt
    x = np.arange(0.5, 6.5, .5)
    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    labels = ['CanESM2', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0',
              'GFDL-CM3', 'GISS-E2-H', 'GISS-E2-R', 'HadGEM2-ES',
              'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'MIROC-ESM', 'Actual']

    ax1 = fig.add_subplot(1, 1, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
    ax1.axhline(j[-1], color='black', linewidth=2, ls='--')

    ax1.scatter(x, j, marker='x', linewidth=3, s=100, color='k',
                label='Full')
    ax1.scatter(x, j_na, marker='x', linewidth=3, s=100, color='r',
                label='Atlantic')
    ax1.scatter(x, j_in, marker='x', linewidth=3, s=100, color='b',
                label='Indian')
    ax1.scatter(x, j_pa, marker='x', linewidth=3, s=100, color='purple',
                label='Pacific')

    ax1.set_ylabel('Jet extension anomaly (ms$^{-1}$)', fontsize=20)
    ax1.set_yticks(np.arange(-1, 4, .5))
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20)
    ax1.set_xticks(x)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    ax1.set_xticklabels(labels, fontsize=20,
                        rotation=60)
    ax1.set_xlim([0.25, 6.25])
    ax1.axvline(5.75, color='k', linewidth=2)
    plt.legend(loc=3, fontsize=16)
    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=-.2, top=.95, bottom=0.2,
                        left=0.1, right=0.9)


def biases(u_c, ncep):
    from scipy import interpolate
    from netcdfread import ncread
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.basemap import Basemap, addcyclic
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lon_ncep = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd.mon.mean.nc', 'lon')
    lat_ncep = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd.mon.mean.nc', 'lat')
    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc', 'lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc', 'lat')

    u_w = u_c[0]
    u_s = u_c[1]
    ncep_s = np.mean(ncep[5:8], axis=0)
    ncep_w = np.mean(ncep[[0, 10, 11]], axis=0)

    h = interpolate.interp2d(lon, lat[::-1], u_w)
    u_w1 = h(lon42, lat42[::-1])
    h = interpolate.interp2d(lon, lat[::-1], u_s)
    u_s1 = h(lon42, lat42[::-1])
    h = interpolate.interp2d(lon_ncep, lat_ncep[::-1], ncep_s)
    ncep_s1 = h(lon42, lat42[::-1])
    h = interpolate.interp2d(lon_ncep, lat_ncep[::-1], ncep_w)
    ncep_w1 = h(lon42, lat42[::-1])

    u_w1, lon1 = addcyclic(u_w1, lon42)
    u_s1, lon1 = addcyclic(u_s1, lon42)
    ncep_s1, lon1 = addcyclic(ncep_s1, lon42)
    ncep_w1, lon1 = addcyclic(ncep_w1, lon42)

    meshlon, meshlat = np.meshgrid(lon1, lat42)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-6, 6, 17)
    ctrs2 = np.concatenate((np.linspace(-20, -2, 10), np.linspace(2, 20, 10)))
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = plt.subplot2grid((3, 6), (0, 0), colspan=5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, u_w1-ncep_w1, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contour(x, y, ncep_w1, ctrs2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(0., 450., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(0., 390., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax2 = plt.subplot2grid((3, 6), (1, 0), colspan=5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=0, urcrnrlon=360, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, u_s1-ncep_s1, ctrs, cmap=newcmap, vmin=np.min(ctrs),
                      vmax=np.max(ctrs), extend='both')
    m.contour(x, y, ncep_s1, ctrs2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    m.drawmeridians(np.arange(0, 450, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(0., 450., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(0., 390., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax4 = plt.subplot2grid((3, 6), (0, 5))
    for axis in ['top', 'bottom', 'left', 'right']:
        ax4.spines[axis].set_linewidth(2)
    ax4.plot(np.mean(u_w1[:, 107:], axis=1), lat42, color='red', label='Model')
    ax4.plot(np.mean(ncep_w1[:, 107:], axis=1), lat42, color='k', label='NCEP2')
    ax4.plot(np.mean(u_w1[:, 107:]-ncep_w1[:, 107:], axis=1), lat42, color='blue', label='Bias')
    ax4.set_ylim([-30, 90])
    ax4.set_xlim([-10, 10])
    ax4.set_yticks(np.arange(0., 90., 30.))
    #ax4.set_xticks(np.arange(-8, 12, 2))
    ax4.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax4.axvline(color='k', ls='--', lw=2)
    ax4.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.9)
    ax4.legend(loc=1)

    ax5 = plt.subplot2grid((3, 6), (1, 5))
    for axis in ['top', 'bottom', 'left', 'right']:
        ax5.spines[axis].set_linewidth(2)
    ax5.plot(np.mean(u_s1[:, 107:], axis=1), lat42, color='red')
    ax5.plot(np.mean(ncep_s1[:, 107:], axis=1), lat42, color='k')
    ax5.plot(np.mean(u_s1[:, 107:]-ncep_s1[:, 107:], axis=1), lat42, color='blue')
    ax5.set_ylim([-30, 90])
    ax5.set_xlim([-10, 10])
    #ax5.set_xticks(np.arange(-8, 12, 2))
    ax5.set_yticks(np.arange(0., 90., 30.))
    ax5.tick_params(labelleft='off', labelbottom='on', which='major',
                    direction='out', length=5, width=2)
    ax5.set_xlabel('Zonal mean u850 (ms$^{-1}$)', size=12)
    ax5.axvline(color='k', ls='--', lw=2)
    ax5.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.9)

    cbar_ax = fig.add_axes([0.233, 0.34, 0.4, 0.005])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='Difference (ms$^{-1}$)', size=20,
                fontsize=12, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    pad = 37.5  # in points

    ax1.annotate('Winter', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax2.annotate('Summer', xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax2.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)

    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=-.2, top=.99, bottom=0.085,
                        left=0.06, right=0.95)


def postage(sst_detrend, gto_interp, nao_ind, month=0, title=''):
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(89.5, -89.5, 180)
    lon = np.linspace(-179.5, 179.5, 360)
    #long = np.concatenate((lon, [180.5]))
    for i in range(7):
        fig, axs = plt.subplots(7, 3, figsize=(15, 14), facecolor='w',
                                edgecolor='k', linewidth=2)
        plotdata0 = sst_detrend[11:-1:12] - sst_detrend[11:-1:12].mean(axis=0)
        plotdata2 = sst_detrend[12::12] - sst_detrend[12::12].mean(axis=0)
        plotdata3 = sst_detrend[13::12] - sst_detrend[13::12].mean(axis=0)
        plotdata = (plotdata0 + plotdata2 + plotdata3)/3
        for a in range(7):
            for b in range(3):
                plotdata1, long = shiftgrid(0., plotdata[(3*a+b)+i*21], lon)
                gto_interp1, long = shiftgrid(0., gto_interp*1e5, lon)
                plotdata1 = np.where(plotdata1==0, np.nan, plotdata1)
                #plotdata1 = np.concatenate((sst_detrend[(a*5+b)*10],
                   #                         np.expand_dims(sst_detrend[(a*5+b)*10, :, 0],
                      #                     axis=1)), axis=1)
                meshlon, meshlat = np.meshgrid(long, lat)
                ax1 = axs[a, b]
                m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                            llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
                m.drawcoastlines()
                m.drawmapboundary(linewidth=2)
                x, y = m(meshlon, meshlat)
                mycmap2 = plt.cm.YlOrRd(np.arange(256))
                mycmap1 = plt.cm.Blues_r(np.arange(256))
                my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
                #my_cmap[239:274, :] = 1
                newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                                       my_cmap)
                caxismin, caxismax, ctrs = -.4,.4, np.linspace(-.4, .4, 17)
                plot = m.contourf(x, y, plotdata1*gto_interp1, ctrs,
                                  cmap=newcmap, vmin=caxismin, vmax=caxismax,
                                  extend='both')
                m.contour(x, y, gto_interp1*(plotdata1*0+1), [-1,-.875,-.75,-.625, -.5, -.375, -.25, -.125, .125,.25,.375,.5,.625,.75,.875,1], colors='k')
                ax1.set_ylim(-60, 60)
                m.drawparallels(np.arange(-60, 90, 30),
                                labels=[False, False, False, False], linewidth=0)
                m.drawmeridians(np.arange(0, 390, 30),
                                labels=[False, False, False, False], linewidth=0)
                ax1.set_yticks(np.arange(-60., 90., 30.))
                ax1.set_xticks(np.arange(0., 390., 30.))
                ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                                direction='out', length=5, width=2)
                ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
                ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
                ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                                direction='out', length=4, width=1)
                ax1.annotate(str((3*a+b)+i*21+1871), xy=(.45, 1.02), fontname='Arial',
                             fontsize=16, xycoords='axes fraction')
                ax1.annotate('{:.2f}'.format(nao_ind[(3*a+b)+i*21]), xy=(.9, 1.02), fontname='Arial',
                             fontsize=16, xycoords='axes fraction')
                '''
                if b == 0:
                    m.drawparallels(np.arange(-60, 90, 30),
                                    labels=[True, False, False, True], linewidth=0)
                if a == 2:
                    m.drawmeridians(np.arange(0, 390, 30),
                                    labels=[True, False, False, True], linewidth=0)
                '''
        cbar_ax = fig.add_axes([0.3, 0.05, 0.4, 0.005])
        b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                         orientation='horizontal', extend='max', format='%.1f')
        b.set_label(label='Sensitivity ([$\mathregular{1x10^5km^2}$]$\mathregular{^{-1}}$)', size=20,
        #b.set_label(label='SST anomaly (K)', size=20,
                    fontsize=12, fontname='Arial', labelpad=-2.5)
        cl = plt.getp(cbar_ax, 'xmajorticklabels')
        plt.setp(cl, fontname='Arial', fontsize=12)
        for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
            label.set_visible(False)
        plt.subplots_adjust(hspace=0.18, wspace=0.1, top=.97, bottom=0.08, left=.01,
                            right=.99)
    
        plt.show()
        plt.savefig('/home/bakerh/Downloads/sstgto_anom_djf_'+ '{:0}'.format(i+1) +'.png',dpi=100)
        plt.close()

def sst_eofs(eof1, eof2, eof3, eof4, winter='yes'):
    from mpl_toolkits.basemap import Basemap, addcyclic
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    eofs = [eof1, eof2, eof3, eof4]
    label = ['a', 'b', 'c', 'd']
    titles = ['EOF1', 'EOF2', 'EOF3', 'EOF4']
    mycmap = plt.cm.RdBu_r(np.arange(256))
    mycmap[115:141, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", mycmap)
    lat = np.linspace(89.5, -89.5, 180)
    lon = np.linspace(-179.5, 179.5, 360)
    if winter == 'yes':
        colormax = .4
        colormin = -colormax
        var = ['47.4', '11.9', '10.4', '6.0'] #detrended
        #var = ['61.0', '9.9', '7.2', '4.1']
    else:
        colormax = .4
        colormin = -colormax
        var = ['44.3', '11.7', '10.0', '6.7'] #  detrended
        var = ['61.8', '7.7', '6.8', '5.2']

    fig, axs = plt.subplots(2, 2, facecolor='w', edgecolor='k', linewidth=2)
    for i in range(4):
        ax1 = axs[int(np.floor(i/2)), np.remainder(i, 2)]
        pdata = eofs[i]
        pdata, lon1 = addcyclic(pdata, lon)
        meshlon, meshlat = np.meshgrid(lon1, lat)
        m = Basemap(projection='cyl', llcrnrlon=-180, llcrnrlat=-60,
                    urcrnrlon=180, urcrnrlat=60, ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        ctrs = np.linspace(colormin, colormax, 17)
        plot = m.contourf(x, y, pdata, ctrs,
                          cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                          extend='both')
        parallels = m.drawparallels(np.arange(-90., 120., 15.))
        meridians = m.drawmeridians(np.arange(-180., 181., 30))
        m.drawparallels(parallels, labels=[True, False, False, False])
        m.drawmeridians(meridians, labels=[False, False, False, True])
        ax1.annotate(label[i], xy=(-.04, .95), fontname='Arial', fontsize=20,
                     fontweight='bold', xycoords='axes fraction')
        ax1.annotate('Var = ' + var[i] + '%', xy=(.78, 1.02), fontname='Arial',
                     fontsize=16, xycoords='axes fraction')
        ax1.set_title(titles[i], fontname='Arial', fontsize=16)

    cbar_ax = fig.add_axes([0.3, 0.08, 0.4, 0.01])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.2f')
    b.set_label(label='EOF (K per SD)', size=16,
                fontsize=16, fontname='Arial', labelpad=0)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    plt.subplots_adjust(hspace=0.2, wspace=0.12, top=.95, bottom=0.15,
                        left=0.03, right=.97)
    plt.show()
